# Kerr Spacetime Geometry and 3+1 ADM Decomposition

## 1. Kerr Metric in Boyer-Lindquist Coordinates

The Kerr metric describes the stationary, axisymmetric vacuum spacetime outside a rotating black hole of mass $M$ and angular momentum $J = a M$. In standard Boyer-Lindquist coordinates $(t, r, \theta, \phi)$, the line element $ds^2 = g_{\mu\nu} dx^\mu dx^\nu$ is given by:

$$ds^2 = -\left(1 - \frac{2Mr}{\Sigma}\right) dt^2 - \frac{4Mar\sin^2\theta}{\Sigma} dt\,d\phi + \frac{\Sigma}{\Delta} dr^2 + \Sigma\,d\theta^2 + \frac{A\sin^2\theta}{\Sigma} d\phi^2$$

where the geometric auxiliary functions are defined as:
$$\Sigma(r, \theta) = r^2 + a^2 \cos^2\theta$$
$$\Delta(r) = r^2 - 2Mr + a^2$$
$$A(r, \theta) = (r^2 + a^2)^2 - a^2 \Delta \sin^2\theta = (r^2 + a^2)\Sigma + 2Mar^2\sin^2\theta$$

### Horizon Radii and Ergosphere
The event horizons correspond to the coordinate singularities where $\Delta(r) = 0$:
$$r_\pm = M \pm \sqrt{M^2 - a^2}$$
where $r_+$ is the outer event horizon and $r_-$ is the inner Cauchy horizon. For physical sub-extremal Kerr black holes, the dimensionless spin parameter satisfies $|a| \le M$.

The stationary limit surfaces (ergosphere boundaries) occur where $g_{tt} = 0$:
$$r_{ergo}(\theta) = M + \sqrt{M^2 - a^2 \cos^2\theta}$$

In the region $r_+ < r < r_{ergo}$, observers cannot remain stationary relative to infinity; they are inevitably dragged in the direction of black hole rotation (frame dragging).

---

## 2. Kerr-Schild Horizon-Penetrating Coordinates

Boyer-Lindquist coordinates exhibit a coordinate divergence at $r = r_+$ because $g_{rr} = \Sigma/\Delta \to \infty$. For numerical simulations of accretion flows falling through the horizon, we employ Kerr-Schild coordinates $(\tilde{t}, x, y, z)$ or spherical Kerr-Schild $(\tilde{t}, r, \theta, \tilde{\phi})$.

The Kerr metric in Kerr-Schild form is written as a Minkowski background perturbed by a null vector field $l_\mu$:
$$g_{\mu\nu} = \eta_{\mu\nu} + 2 H l_\mu l_\nu$$
where $\eta_{\mu\nu} = \text{diag}(-1, 1, 1, 1)$ is the flat Minkowski metric, and the scalar function is:
$$H(r, \theta) = \frac{M r^3}{r^4 + a^2 z^2}$$
The vector $l^\mu$ is null and geodesic with respect to both $\eta_{\mu\nu}$ and $g_{\mu\nu}$:
$$g^{\mu\nu} l_\mu l_\nu = \eta^{\mu\nu} l_\mu l_\nu = 0$$

Because $\det(g_{\mu\nu}) = \det(\eta_{\mu\nu}) = -1$ in Cartesian Kerr-Schild coordinates, the metric determinant remains finite, smooth, and well-behaved across the horizon $r = r_+$.

---

## 3. 3+1 ADM Foliation for General Relativistic MHD

To formulate hyperbolic evolution equations, spacetime $(\mathcal{M}, g)$ is foliated into a family of spacelike hypersurfaces $\Sigma_t$ parameterized by coordinate time $t$.

The 4-metric line element decomposes into the ADM form:
$$ds^2 = -\alpha^2 dt^2 + \gamma_{ij} (dx^i + \beta^i dt)(dx^j + \beta^j dt)$$

where:
- $\alpha$: Lapse function, measuring the rate of advance of proper time relative to coordinate time along hypersurface normals.
- $\beta^i$: Shift vector, measuring the tangential dragging of spatial coordinates.
- $\gamma_{ij}$: Induced 3-metric on the spatial hypersurface $\Sigma_t$.

The inverse metric components are:
$$g^{00} = -\frac{1}{\alpha^2}$$
$$g^{0i} = \frac{\beta^i}{\alpha^2}$$
$$g^{ij} = \gamma^{ij} - \frac{\beta^i \beta^j}{\alpha^2}$$

The 4-metric determinant relates to the 3-metric determinant via:
$$\sqrt{-g} = \alpha \sqrt{\gamma}$$

---

## 4. Null Geodesics and the Bardeen Shadow Boundary

Photons follow null geodesics $ds^2 = g_{\mu\nu} \dot{x}^\mu \dot{x}^\nu = 0$. In Kerr spacetime, Carter (1968) demonstrated that the Hamilton-Jacobi equation separates completely due to the presence of a second-order Killing tensor $K_{\mu\nu}$:
$$\nabla_{(\lambda} K_{\mu\nu)} = 0$$

The four constants of motion are:
1. $E = -p_t$: Energy at infinity.
2. $L_z = p_\phi$: Axial angular momentum.
3. $\mu^2 = -p_\mu p^\mu = 0$: Rest mass.
4. $Q = p_\theta^2 + \cos^2\theta \left( a^2 E^2 + \frac{L_z^2}{\sin^2\theta} \right)$: Carter constant.

Introducing impact parameters $\xi = L_z / E$ and $\eta = Q / E^2$, the radial equation of motion reads:
$$\Sigma \frac{dr}{d\lambda} = \pm \sqrt{R(r)}$$
$$R(r) = (r^2 + a^2 - a\xi)^2 - \Delta \left( \eta + (\xi - a)^2 \right)$$

Unstable circular photon orbits satisfying $R(r_c) = 0$ and $R'(r_c) = 0$ define the boundary of the black hole shadow. Solving these simultaneous relations yields the analytical Bardeen parameterization for the critical impact parameters:
$$\xi_c(r_c) = \frac{(r_c^2 + a^2)\Delta'(r_c) - 4 r_c \Delta(r_c)}{a \Delta'(r_c)}$$
$$\eta_c(r_c) = \frac{16 r_c^2 \Delta(r_c) - \left( (r_c^2 - a^2)\Delta'(r_c) - 4 r_c \Delta(r_c) \right)^2}{a^2 [\Delta'(r_c)]^2}$$

For an observer at inclination $\theta_o$ and large distance $r_o \gg M$, the celestial coordinates $(\alpha_{cel}, \beta_{cel})$ on the image plane map to:
$$\alpha_{cel} = -\frac{\xi_c}{\sin\theta_o}$$
$$\beta_{cel} = \pm \sqrt{\eta_c + a^2 \cos^2\theta_o - \xi_c^2 \cot^2\theta_o}$$

In the Schwarzschild limit ($a \to 0$), the shadow is a perfect circle of radius $R_{shadow} = \sqrt{27} M \approx 5.196152 M$. For spinning black holes, the shadow is deformed asymmetrically and displaced along the horizontal axis due to frame dragging.
