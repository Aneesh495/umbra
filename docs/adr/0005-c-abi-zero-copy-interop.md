# ADR 0005: Zero-Copy C ABI and Python Foreign Function Interface

## Status
Accepted

## Context
Scientific computational engines written in Fortran must interface seamlessly with modern data science, visualization, and machine learning ecosystems (such as NumPy, PyTorch, SciPy, and Matplotlib). Traditional interfacing approaches (such as writing intermediate disk files or using legacy F77 wrappers) introduce severe I/O bottlenecks and memory copies.

To ensure frictionless adoption, `umbra` required a zero-copy foreign function interface (FFI) adhering strictly to standard ISO C interoperability without external tool dependencies (such as f2py or SWIG).

## Decision
We establish a native C Application Binary Interface (ABI) directly in Fortran 2003/2018 using `iso_c_binding`:

1. **Standard ABI Declarations**:
   All public API functions are declared with `bind(C, name="...")` and utilize standardized C types (`c_int`, `c_double`, `c_ptr`, `c_char`).
2. **Opaque Pointer Encapsulation**:
   Complex Fortran derived types (such as `mesh_t`, `dg_mesh_t`, `geodesic_camera_t`) are encapsulated behind opaque C pointer handles (`c_ptr`). Fortran's internal type state is managed on the heap using pointers and converted via `c_f_pointer`.
3. **Contiguous Zero-Copy Memory Pointers**:
   State vectors and image buffers are exposed directly as flat contiguous double-precision arrays (`real(c_double), dimension(*)`) without intermediate memory allocation or copying.
4. **Header Architecture**:
   A single clean C header file `include/umbra.h` provides complete ANSI C function signatures, struct layouts, and return error codes.
5. **Python ctypes Binding Layer**:
   A pure-Python wrapper package `bindings/python/umbra/` binds directly to `libumbra.so` / `libumbra.dylib` using standard library `ctypes`. Python `numpy.ctypeslib` provides zero-copy buffer views over native Fortran mesh allocations.

## Consequences
- **Advantages**:
  - Zero computational or memory overhead for data passing between Fortran and external host languages (C, C++, Python, Julia, Rust).
  - Eliminates build-time dependencies on wrapper generators; compiles cleanly with any standard C/Fortran compiler.
  - Python users gain interactive REPL workflows, Matplotlib visualizations, and automated test orchestration while preserving raw Fortran kernel performance.
- **Trade-offs**:
  - Memory lifecycle management requires explicit allocation and deallocation helper routines (`umbra_mesh_free`, `umbra_camera_free`) to prevent memory leaks across language boundaries.
