# ADR 0001: Nodal Discontinuous Galerkin Polynomial Basis

## Status
Accepted

## Context
General Relativistic Magnetohydrodynamics (GRMHD) in curved spacetimes exhibits steep gradients, relativistic shocks, and turbulent magnetorotational instabilities (MRI). Classical high-order finite difference and finite volume schemes require wide stencils across multiple mesh cells, introducing boundary communication bottlenecks and stencil truncation near coordinate singularities (such as coordinate cutoffs or black hole event horizons).

To achieve high-order accuracy with compact computational stencils, we required a spatial discretization scheme that supports element-local high-order representation and natural mesh refinement.

## Decision
We adopt a Nodal Discontinuous Galerkin (DG) method using Legendre-Gauss-Lobatto (LGL) quadrature nodes and Lagrange interpolating polynomials:

1. Spatial coordinates within each reference element $[-1, 1]^d$ are mapped to physical coordinates using affine or curvilinear geometric Jacobians.
2. Solution states are represented nodally at the tensor-product LGL points.
3. Integration of volume terms and boundary fluxes uses the same LGL quadrature rule. Under LGL quadrature, the Lagrange polynomials are mutually orthogonal, resulting in a diagonal mass matrix (exact collocation property).
4. Boundary inter-element fluxes are mediated by approximate Riemann solvers (such as HLL, HLLC, or HLLD) evaluated at collocated surface quadrature nodes.
5. High-order oscillations near shock discontinuities are controlled via positivity-preserving limiters (such as Zhang-Shu limiter) and minmod slope limiters.

## Consequences
- **Advantages**:
  - The mass matrix is strictly diagonal, eliminating the need to invert dense global or element-level mass matrices at each Runge-Kutta sub-stage.
  - Computational stencils only require immediate nearest-neighbor element face communication, enabling ideal strong scaling on distributed memory architectures.
  - Sub-element resolution enables accurate capture of turbulent eddy cascades without excessive numerical dissipation.
- **Trade-offs**:
  - Severe Courant-Friedrichs-Lewy (CFL) time step restriction proportional to $1 / (2k + 1)$ for polynomial order $k$.
  - Strong shocks can trigger Gibbs phenomena unless bounded by robust limiters and physical entropy floors.
