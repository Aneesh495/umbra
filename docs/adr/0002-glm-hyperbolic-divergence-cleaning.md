# ADR 0002: Hyperbolic Divergence Cleaning via Generalized Lagrange Multiplier (GLM)

## Status
Accepted

## Context
In classical and relativistic magnetohydrodynamics, Gauss's law for magnetism mandates the solenoidal constraint:

$$\nabla \cdot \mathbf{B} = 0$$

In curved spacetime 3+1 formulations, this condition translates to $\nabla_i B^i = 0$ on the spatial hypersurface $\Sigma_t$. Standard collocated shock-capturing schemes do not automatically preserve this constraint. Numerical accumulation of non-zero magnetic monopoles triggers spurious magnetic forces along field lines, corrupts the relativistic Lorentz factor, and causes catastrophic code crashes near event horizons.

Two primary strategies exist in the literature:
1. Constrained Transport (CT) on staggered face-centered meshes (Evans & Hawley, Stone & Gardiner).
2. Hyperbolic Divergence Cleaning via a scalar potential $\psi$ (Dedner et al., 2002).

## Decision
We implement a dual approach with primary support for cell-centered Hyperbolic Divergence Cleaning (Dedner GLM) augmented with staggered Constrained Transport:

1. The magnetic induction equation is coupled to an auxiliary scalar field $\psi$:
   $$\partial_t B^i + \dots + \partial_i \psi = 0$$
   $$\partial_t \psi + c_h^2 \partial_i B^i = -\frac{c_h^2}{c_p^2} \psi$$
2. Monopole errors are propagated hyperbolically across the domain at speed $c_h$ (equal to the maximum characteristic speed, typically $c=1$) toward boundary sinks.
3. Monopole errors are simultaneously damped exponentially at rate $\tau = c_p^2 / c_h^2$.
4. For curvilinear Kerr-Schild metrics, the coupling is modified to incorporate the spatial metric determinant $\gamma$ and lapse function $\alpha$.

## Consequences
- **Advantages**:
  - Operates directly on cell-centered or nodal DG state vectors without requiring staggered meshes or staggered edge interpolations.
  - Seamlessly integrates with standard Riemann solvers (such as HLLD and HLL) and high-order Runge-Kutta time integrators.
  - Enables local adaptive mesh refinement (AMR) without complex hanging-edge divergence correction loops.
- **Trade-offs**:
  - Expands the conservative state vector from 8 variables $(\rho^*, \tilde{S}_i, \tau, \tilde{B}^i)$ to 9 variables with $\psi$.
  - Damping parameter $c_p$ must be calibrated against the grid scale to balance wave propagation against over-damping.
