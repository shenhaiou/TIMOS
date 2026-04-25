# TIM-OS (Tetrahedron-based Inhomogeneous Monte-Carlo Optical Simulator)

**TIM-OS** is a high-performance Monte Carlo photon transport simulator designed for biomedical optics. Originally developed by **Haiou Shen** and **Ge Wang** (2009-2010), it provides a flexible and efficient framework for modeling light propagation in complex, inhomogeneous biological systems using unstructured tetrahedral meshes.

### Key Features
- **Complex Geometries:** Uses tetrahedral meshes to accurately represent irregular boundaries and anatomical structures (e.g., small animal models or human organs).
- **Inhomogeneity Support:** Assigns independent optical properties (absorption, scattering, anisotropy, and refractive index) to each tetrahedral element.
- **Modern C++23 Core:** Completely modernized (v3.5) with strict RAII, `std::vector` memory management, and `std::expected` error handling.
- **High Performance:**
  - **SIMD Optimized:** Utilizes ARM NEON (Apple Silicon) and scalar fallbacks for fast ray-tetrahedron intersection tests.
  - **Hardware Acceleration:** Leverages the Apple Accelerate framework for cryptographically secure, high-throughput AES-CTR random number generation.
  - **Scalable Multi-threading:** Uses `std::jthread` for efficient cross-platform parallel execution.
- **Benchmark Mode:** Includes a pure compute mode for measuring CPU and memory performance without I/O bottlenecks.

### Scientific Foundation
This simulator is based on the methods detailed in the following original publications:
1.  **H. Shen and G. Wang**, "A tetrahedron-based inhomogeneous Monte Carlo optical simulator," *Physics in Medicine and Biology*, vol. 55, no. 4, pp. 947-962, 2010.
2.  **H. Shen and G. Wang**, "A study on tetrahedron-based inhomogeneous Monte Carlo optical simulation," *Biomedical Optics Express*, vol. 1, no. 1, pp. 110-123, 2010.

## Build and Usage

### Requirements
- A modern C++ compiler supporting **C++23** (e.g., Clang 16+, GCC 13+).
- (Optional) Apple Silicon Mac with the Accelerate framework for maximum performance.

### Build
On macOS/Apple Silicon:
```bash
make       # Release build with LTO
make pgo   # Profile-guided optimization build (~20% faster)
```

### Run
```bash
./source/timos -f mesh_file -s source_file -p optics_file -m [s|i|si] -o output_file [options]
```
For benchmarking without disk I/O, simply omit the `-o` flag.

## License and Attribution
- **Original Development (2009):** Haiou Shen (Virginia Tech)
- **Modernization (2026):** Haiou Shen
