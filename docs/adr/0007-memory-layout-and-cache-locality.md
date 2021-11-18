# ADR 0007: Column-Major Memory Layout and Cache Locality

## Status
Accepted

## Context
High-performance PDE kernels and radiative transfer engines are frequently memory-bandwidth bound rather than arithmetic-logic bound. In Fortran, multidimensional arrays are stored in column-major order (the first index varies fastest in contiguous memory addresses).

Non-contiguous striding (traversing the rightmost index in the innermost loop) causes cache line thrashing, instruction cache misses, and disables SIMD auto-vectorization. Furthermore, discontinuous Galerkin methods evaluate small tensor contractions at each quadrature node, demanding spatial and temporal locality.

## Decision
We architect all memory layouts to align with Fortran's column-major ordering and L1/L2 cache line geometries:

1. **Inner Loop Indexing**:
   All nested loops iterate with the first index in the innermost loop:
   ```fortran
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            u(i, j, k, var) = ...
         end do
      end do
   end do
   ```
2. **Contiguous Array Slicing**:
   When extracting state vectors at element boundaries or cell interfaces, contiguous slices `u(:, j, k)` are utilized to ensure compiler vectorization directives (`!$omp simd`) can generate packed SIMD instructions.
3. **Structure of Arrays (SoA) vs Array of Structures (AoS)**:
   For bulk grid storage, multidimensional state arrays `u(nx, ny, nz, nvar)` place spatial indices first and physical variable indices last, allowing vectorized operations across adjacent spatial cells for a given physical component.
4. **Local Quadrature Scratch Buffers**:
   Within Discontinuous Galerkin volume and surface routines, nodal values are mapped into small stack-allocated or thread-local contiguous buffers that fit entirely within L1 data cache (typically 32 KB to 64 KB per core).

## Consequences
- **Advantages**:
  - Maximum cache line utilization (100% cache line hit rate on linear sweeps).
  - Clean automatic vectorization by compilers (GFortran `-O3`, LLVM Flang, Intel ifx).
  - Eliminates unneeded cache evictions during DG matrix-vector contractions.
- **Trade-offs**:
  - Requires developers to maintain discipline regarding loop ordering conventions across multi-dimensional tensors.
