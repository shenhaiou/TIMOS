# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**TIM-OS (Tetrahedral-mesh based Inhomogeneous Monte-Carlo Optical Simulator)** v3.0 — a Monte Carlo photon transport simulator for biomedical optics. It propagates photons through a 3D finite element mesh of tetrahedral elements, each with independently assigned optical properties, and records surface and/or internal fluence.

## Build

On **Apple Silicon (M1/M2/M3)** — use the Makefile:

```sh
make          # fast release build (development)
make pgo      # profile-guided build — ~20% faster, takes ~3min
```

The compiled binary is `source/timos`.

**What the optimised build uses:**
- `-march=native -ffast-math -flto` — M1-specific codegen + LTO
- `-DHAS_ACCELERATE -framework Accelerate` — `BNNSRandomGenerator` (AES-CTR) + `vvlog`/`vvsincos` for batched RNG math
- ARM NEON (`arm_neon.h`) for `PhotonTetrahedronIntersection` — all 4 dot products in 2 NEON FMA pairs with sequential `vld1q_f64` from SoA `TriNorm` layout
- `os_unfair_lock` instead of `pthread_mutex` — lower-overhead Apple Silicon lock
- Per-thread local accumulators (replace the 12 MB `PhotonInfo` buffer) — one final lock per thread instead of per-batch locking

On **Intel** (original): `icpc -static -mkl -O3 -o source/timos source/TIMOS.cpp`

## Running a Simulation

```sh
./source/timos \
  -f <mesh_file>    \   # finite element mesh (.mesh)
  -s <source_file>  \   # source definition (.source)
  -p <optics_file>  \   # optical parameters (.opt)
  -m <mode>         \   # s=surface, i=internal, is/si=both
  -o <output_file>  \   # output filename (.dat)
  -t <num_threads>      # thread count (use num_cores or 2x num_cores)
```

Optional flags:
- `-T <time_step_ns> <num_steps>` — enable time-domain simulation
- `-r <start_rand_id>` — starting random stream index (for distributed runs across machines)
- `-k <mask_file>` — mask file

Run all provided examples (32 threads each):
```sh
sh runall.sh
```

## Input File Formats

**`.mesh`** — finite element mesh:
```
numNodes
numElems
x y z          # one line per node (1-indexed)
...
n1 n2 n3 n4 medIdx   # one line per tetrahedron (nodes 1-indexed, medIdx 1-indexed)
...
```

**`.opt`** — optical parameters:
```
1                          # always 1 (file version flag)
numMedia
mua mus g refractiveIndex  # one line per medium
1                          # boundary flag: 0=uniform, 1=matched
envRefractiveIndex
```

**`.source`** — light source(s):
```
numSources
sourceType elemIdx [params...] numPhotons
```
Source types: `1`=isotropic point, `2`=isotropic region, `11`=pencil beam, `12`=surface triangle region.

## Code Architecture

The codebase is a single-file C++ implementation (`source/TIMOS.cpp`, ~2400 lines) plus a header (`source/timos.h`).

**Key data structures** (defined in `timos.h`):
- `TNode` — 3D coordinate/vector
- `TElem` / `TElemNode` — tetrahedral element: stores medium index, four adjacent element indices (negative = boundary triangle), and inward-pointing triangle normal vectors
- `TTriangle` / `TTriNode` — triangle with plane equation (`ax+by+cz+d=0`) and adjacency to 1–2 elements
- `TMedOptic` — optical properties: `mua`, `mus`, `g`, refractive index, plus precomputed Henyey-Greenstein terms
- `TSource` — source configuration
- `TPhoton` — per-photon state: position, direction, weight, current element

**Global state** (all prefixed `g_`):
- `g_Nodes`, `g_Elems`, `g_ElemNodes`, `g_Triangles`, `g_BoundaryTrigs` — mesh
- `g_MedOptic` — optical parameters per medium (or per element if `g_SimpleOptic=false`)
- `g_SurfMeas`, `g_Absorption` — accumulated fluence results
- `g_TimeDomain`, `g_TimeStep`, `g_NumTimeStep` — time-domain settings

**Threading**: pthreads with `Result_Lock` and `Source_Lock` mutexes. Each thread uses an independent Intel MT2203 random stream (up to 1024 streams), enabling distributed computation across machines by setting `-r` to a non-overlapping starting index.

**Simulation flow**:
1. Parse arguments → load mesh (`fem_read`) → build triangle adjacency structure
2. Load optical parameters (`optic_read`) → precompute Henyey-Greenstein coefficients
3. Load sources (`source_read`) → locate source element in mesh
4. Launch `g_NumThread` pthreads, each running the photon loop
5. Each photon: sample step length, walk through tetrahedra via adjacency, handle boundary reflection/refraction (Fresnel), apply absorption weight
6. Accumulate results under mutex → write output

**Distributed/cluster use**: split total photons across machines with different `-r` values. Each machine's thread count determines the stride (e.g., 16 threads/machine → machine 1: `-r 1`, machine 2: `-r 17`).
