# TIM-OS (Tetrahedral-mesh based Inhomogeneous Monte-Carlo Optical Simulator)

TIM-OS v3.0 is a high-performance Monte Carlo photon transport simulator designed for biomedical optics. It propagates photons through 3D finite element meshes composed of tetrahedral elements, allowing for complex inhomogeneous geometries.

## Project Overview

- **Purpose:** Simulate light propagation in biological tissues using tetrahedral meshes.
- **Main Technologies:** C++23, ARM NEON SIMD (Apple Silicon), Apple Accelerate Framework, PGO (Profile-Guided Optimization).
- **Architecture:** Refactored into modular components (Mesh, Optics, Source, Simulation, Output) managed by a central `SimContext`.

## Key Components

- **`source/mesh.cpp` / `mesh.h`**: Handles loading `.mesh` files and building the adjacency structure for tetrahedra and triangles.
- **`source/optics.cpp` / `optics.h`**: Manages optical properties (mua, mus, g, refractive index) and precomputes scattering parameters.
- **`source/source_io.cpp` / `source_io.h`**: Loads light source definitions (point, beam, etc.) and initializes photons.
- **`source/simulation.cpp` / `simulation.h`**: The core Monte Carlo engine, implementing photon stepping, boundary interaction (Fresnel), and absorption.
- **`source/output.cpp` / `output.h`**: Records fluence results and writes them to disk in ASCII format.
- **`source/simcontext.h`**: Central state container using RAII (std::vector) to manage mesh and result data.

## Building and Running

### Build Commands
- `make`: Standard release build with LTO (Link Time Optimization).
- `make pgo`: Optimized build using Profile-Guided Optimization (~20% faster on Apple Silicon).
- `make test`: Compiles and runs the test suite.
- `make clean`: Removes binaries and profile artifacts.

### Running a Simulation
```bash
./source/timos -p <opt_file> -f <mesh_file> -s <source_file> -m <mode> -o <output_file> [options]
```
- **Required Flags:**
  - `-p, --optical`: Optical parameter file (`.opt`).
  - `-f, --mesh`: Finite element mesh file (`.mesh`).
  - `-s, --source`: Light source definition (`.source`).
  - `-m, --mode`: Output mode (`s`=surface, `i`=internal, `si`=both).
  - `-o, --output`: Output data filename.
- **Optional Flags:**
  - `-t, --threads N`: Number of worker threads (default: 1).
  - `-r, --rand-start N`: Starting RNG stream index (useful for distributed cluster runs).
  - `-T, --time-step DT`: Enable time-domain mode with step size DT in ns.
  - `-N, --num-steps N`: Number of time steps for time-domain mode.

## Development Conventions

- **C++ Standards:** Uses C++23 features (e.g., `std::expected`, `std::jthread`).
- **Memory Management:** Prefers RAII and `std::vector` over manual pointers. Arrays are typically 1-indexed to match legacy mesh formats (index 0 is unused).
- **Performance:** Leverages ARM NEON intrinsics for intersection math and the Accelerate framework for fast RNG (AES-CTR).
- **Testing:** Integration and unit tests are located in the `tests/` directory. New features should include corresponding test cases in `test_integration.cpp` or specialized test files.
- **Concurrency:** Uses `std::jthread` for thread management and atomic/mutex-protected accumulation for results.

## Input/Output Formats

Refer to `CLAUDE.md` for detailed specifications of `.mesh`, `.opt`, and `.source` file formats. Output is typically saved as `.dat` files.
