/*
 * umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
 * Header: umbra.h
 * Description: C API bindings and ABI declarations for foreign runtime interop.
 */

#ifndef UMBRA_H
#define UMBRA_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Library version info */
int32_t umbra_version_major(void);
int32_t umbra_version_minor(void);
int32_t umbra_version_patch(void);

/* Spacetime & Shadow Geometry */
double umbra_kerr_shadow_schwarzschild_radius(void);

void umbra_kerr_shadow_analytic(
    double spin,
    double theta_obs,
    int32_t n_points,
    double* alpha_out,
    double* beta_plus_out,
    double* beta_minus_out
);

/* Conservative-to-Primitive Inversion */
int32_t umbra_c_con2prim_1d(
    double d, double s1, double s2, double s3, double tau,
    double b1, double b2, double b3,
    double gamma_eos,
    double* rho_out, double* vx_out, double* vy_out, double* vz_out, double* p_out
);

int32_t umbra_c_con2prim_2d(
    double d, double s1, double s2, double s3, double tau,
    double b1, double b2, double b3,
    double gamma_eos,
    double* rho_out, double* vx_out, double* vy_out, double* vz_out, double* p_out
);

/* Relativistic Riemann Solvers */
void umbra_c_riemann_llf(
    const double* prim_l,
    const double* prim_r,
    double gamma_eos,
    int32_t dir,
    double* flux_out
);

void umbra_c_riemann_hll(
    const double* prim_l,
    const double* prim_r,
    double gamma_eos,
    int32_t dir,
    double* flux_out
);

void umbra_c_riemann_hlld(
    const double* prim_l,
    const double* prim_r,
    double gamma_eos,
    int32_t dir,
    double* flux_out
);

/* Null Geodesic Step */
int32_t umbra_c_geodesic_step(
    double spin,
    double d_lambda,
    double r_horizon,
    double r_escape,
    double* x,
    double* p,
    bool* terminated,
    bool* captured
);

/* Camera Raytracing */
void umbra_c_render_shadow(
    double spin,
    double inclination_deg,
    double fov_m,
    int32_t npix_x,
    int32_t npix_y,
    int32_t* mask_out
);

#ifdef __cplusplus
}
#endif

#endif /* UMBRA_H */
