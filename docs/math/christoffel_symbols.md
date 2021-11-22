# Christoffel Symbols and Connection Geometries

## 1. Metric Connection Definition

On a pseudo-Riemannian spacetime manifold equipped with a metric tensor $g_{\mu\nu}$, the unique symmetric, metric-compatible Levi-Civita connection coefficients (Christoffel symbols of the second kind) are defined as:

$$\Gamma^\lambda_{\mu\nu} = \frac{1}{2} g^{\lambda\sigma} \left( \partial_\mu g_{\nu\sigma} + \partial_\nu g_{\mu\sigma} - \partial_\sigma g_{\mu\nu} \right)$$

where $\partial_\sigma = \frac{\partial}{\partial x^\sigma}$.

The symmetry property holds identically under index exchange:
$$\Gamma^\lambda_{\mu\nu} = \Gamma^\lambda_{\nu\mu}$$

---

## 2. Geodesic Trajectory Formulation

The covariant derivative of a 4-velocity $u^\mu = dx^\mu / d\tau$ along its own path defines the geodesic equation:

$$\nabla_u u^\lambda = \frac{d^2 x^\lambda}{d\tau^2} + \Gamma^\lambda_{\mu\nu} \frac{dx^\mu}{d\tau} \frac{dx^\nu}{d\tau} = 0$$

For null geodesics parameterized by an affine parameter $\lambda$, this equation describes photon paths in curved spacetime. In numerical ray-tracing, solving this second-order system directly is mathematically equivalent to the first-order Hamiltonian system:

$$\frac{dx^\mu}{d\lambda} = g^{\mu\nu} p_\nu$$
$$\frac{dp_\mu}{d\lambda} = -\frac{1}{2} \left( \partial_\mu g^{\alpha\beta} \right) p_\alpha p_\beta$$

Notice that metric inverse derivatives relate to metric derivatives via:
$$\partial_\mu g^{\alpha\beta} = -g^{\alpha\rho} g^{\beta\sigma} \partial_\mu g_{\rho\sigma}$$

Using the identity $p_\lambda = g_{\lambda\nu} u^\nu$, the equivalence is established directly:
$$\frac{dp_\lambda}{d\lambda} = g_{\lambda\nu} \frac{du^\nu}{d\lambda} + u^\nu \frac{\partial g_{\lambda\nu}}{\partial x^\sigma} u^\sigma = \Gamma^\nu_{\rho\sigma} g_{\lambda\nu} u^\rho u^\sigma = \frac{1}{2} (\partial_\lambda g_{\mu\nu}) u^\mu u^\nu$$

---

## 3. Christoffel Symbols in 3+1 ADM Evolution

In the 3+1 formulation of General Relativistic Hydrodynamics, spatial covariant derivatives $\mathcal{D}_i$ on the hypersurface $\Sigma_t$ are evaluated using the 3-Christoffel connection:

$$^{(3)}\Gamma^k_{ij} = \frac{1}{2} \gamma^{kl} \left( \partial_i \gamma_{jl} + \partial_j \gamma_{il} - \partial_l \gamma_{ij} \right)$$

The geometric source terms in the conservation of momentum equation:
$$\partial_t \tilde{S}_j + \partial_i \mathbf{F}^i_j = \mathbf{S}_j$$
are driven by spatial gradients of the lapse function $\partial_j \alpha$, the shift vector $\partial_j \beta^i$, and the metric tensor:

$$\mathbf{S}_j = \frac{1}{2} \alpha \sqrt{\gamma} T^{\mu\nu} \partial_j g_{\mu\nu}$$

where $T^{\mu\nu} = (w + b^2) u^\mu u^\nu + \left(p + \frac{1}{2} b^2\right) g^{\mu\nu} - b^\mu b^\nu$ is the ideal relativistic MHD stress-energy tensor.

---

## 4. Analytical Derivatives in Boyer-Lindquist Coordinates

For Boyer-Lindquist coordinates, non-vanishing derivatives with respect to $r$ and $\theta$ can be evaluated in closed analytic form:

$$\frac{\partial \Sigma}{\partial r} = 2r, \quad \frac{\partial \Sigma}{\partial \theta} = -2 a^2 \sin\theta \cos\theta$$
$$\frac{\partial \Delta}{\partial r} = 2r - 2M, \quad \frac{\partial \Delta}{\partial \theta} = 0$$

$$\frac{\partial g_{tt}}{\partial r} = \frac{2M(r^2 - a^2\cos^2\theta)}{\Sigma^2}$$
$$\frac{\partial g_{tt}}{\partial \theta} = -\frac{4 M a^2 r \sin\theta \cos\theta}{\Sigma^2}$$

$$\frac{\partial g_{t\phi}}{\partial r} = -\frac{2 M a \sin^2\theta (a^2\cos^2\theta - r^2)}{\Sigma^2}$$
$$\frac{\partial g_{t\phi}}{\partial \theta} = -\frac{4 M a r \sin\theta \cos\theta (r^2 + a^2)}{\Sigma^2}$$

$$\frac{\partial g_{rr}}{\partial r} = \frac{2r}{\Delta} - \frac{\Sigma (2r - 2M)}{\Delta^2}$$
$$\frac{\partial g_{rr}}{\partial \theta} = -\frac{2 a^2 \sin\theta \cos\theta}{\Delta}$$

Evaluating these derivatives analytically rather than by numerical finite differences guarantees machine-precision conservation of the Hamiltonian constraint across millions of numerical integration steps.
