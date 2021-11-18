# ADR 0006: IEEE-754 64-Bit Double Precision Representation

## Status
Accepted

## Context
General Relativistic Magnetohydrodynamics and geodesic ray-tracing span extraordinary dynamic ranges:
- Spacetime curvature gradients $\sim M/r^3$ near black hole horizons.
- Ultra-relativistic Lorentz factors $\Gamma \sim 50$ where $1 - v \sim 2 \times 10^{-4}$.
- Magnetization ratios $\sigma = B^2 / \rho \gg 10^3$.
- Ray trajectories accumulating $10^5$ integration steps with invariant conservation thresholds $\le 10^{-6}$.

Using single-precision 32-bit floating-point numbers (`real32`) causes severe catastrophic cancellation in metric determinants $\sqrt{-g} = \alpha \sqrt{\gamma}$, Lorentz contraction evaluations $\Gamma = 1/\sqrt{1 - v^2}$, and Carter invariant drift. Conversely, quadruple precision (`real128`) incurs a 10x to 20x performance penalty on modern SIMD hardware without hardware support.

## Decision
We enforce double-precision IEEE-754 floating-point arithmetic throughout the entire codebase:

1. **Explicit Kind Specification**:
   All real numbers are declared using the standard Fortran 2008 intrinsic module `iso_fortran_env: real64` aliased to `dp`:
   ```fortran
   use, intrinsic :: iso_fortran_env, only: real64, int32, int64
   integer, parameter :: dp = real64
   ```
2. **Explicit Literal Constants**:
   Every floating-point constant in mathematical expressions includes the explicit `_dp` suffix (for example: `1.0_dp / 3.0_dp`, `0.5_dp`, `1.0e-12_dp`). Bare literals (such as `1.0` or `0.5`) are disallowed to prevent implicit single-precision conversion during compiler intermediate representation generation.
3. **Machine Precision Tolerances**:
   Convergence thresholds and floor values are established relative to double-precision machine epsilon:
   $$\epsilon_{DP} = 2^{-52} \approx 2.2204 \times 10^{-16}$$
   - con2prim Newton-Raphson tolerance: $10^{-12}$
   - Hamiltonian geodesic invariant tolerance: $10^{-6}$
   - Carter constant conservation tolerance: $10^{-4}$

## Consequences
- **Advantages**:
  - Consistent bitwise reproducibility across diverse architectures (x86_64, AArch64/ARM64 Apple Silicon) adhering to IEEE-754 standards.
  - Optimal hardware vectorization on modern AVX-512 and ARM NEON SIMD instruction pipelines.
  - High dynamic range preventing catastrophic cancellation near ergosphere horizons.
- **Trade-offs**:
  - Higher memory bandwidth utilization compared to single precision; mitigated through cache-friendly struct-of-arrays memory layouts.
