! umbra: High-Order Discontinuous Galerkin GRMHD Kernel & Kerr Radiative Transfer
! Module: umbra_math
! Description: Linear algebra, tensor operations, 3x3/4x4 inversions, and nonlinear root solvers.

module umbra_math
  use umbra_kinds, only: dp, i4
  use umbra_constants, only: EPSILON_DP, TOL_CONVERGENCE, PI, TWO_PI
  use umbra_error, only: umbra_error_t, UMBRA_SUCCESS, UMBRA_ERR_SINGULAR_METRIC, &
                         UMBRA_ERR_CON2PRIM_FAILED
  implicit none
  private

  public :: vec3_dot, vec3_cross, vec3_norm, vec3_normalize
  public :: mat3_det, mat3_inv, mat3_vec_mul, mat3_transpose
  public :: mat4_det, mat4_inv, mat4_vec_mul
  public :: mat3_symm_eigenvalues
  public :: solve_newton_bisection
  public :: solve_brent

  abstract interface
    function scalar_func_1d(x) result(f)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: f
    end function scalar_func_1d

    subroutine scalar_func_and_deriv_1d(x, f, df)
      import :: dp
      real(dp), intent(in) :: x
      real(dp), intent(out) :: f
      real(dp), intent(out) :: df
    end subroutine scalar_func_and_deriv_1d
  end interface

contains

  pure function vec3_dot(a, b) result(res)
    real(dp), dimension(3), intent(in) :: a, b
    real(dp) :: res
    res = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end function vec3_dot

  pure function vec3_cross(a, b) result(res)
    real(dp), dimension(3), intent(in) :: a, b
    real(dp), dimension(3) :: res
    res(1) = a(2)*b(3) - a(3)*b(2)
    res(2) = a(3)*b(1) - a(1)*b(3)
    res(3) = a(1)*b(2) - a(2)*b(1)
  end function vec3_cross

  pure function vec3_norm(a) result(res)
    real(dp), dimension(3), intent(in) :: a
    real(dp) :: res
    res = sqrt(vec3_dot(a, a))
  end function vec3_norm

  pure function vec3_normalize(a) result(res)
    real(dp), dimension(3), intent(in) :: a
    real(dp), dimension(3) :: res
    real(dp) :: n
    n = vec3_norm(a)
    if (n > 1.0e-30_dp) then
      res = a / n
    else
      res = 0.0_dp
    end if
  end function vec3_normalize

  pure function mat3_det(a) result(det)
    real(dp), dimension(3,3), intent(in) :: a
    real(dp) :: det
    det = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
  end function mat3_det

  pure function mat3_transpose(a) result(res)
    real(dp), dimension(3,3), intent(in) :: a
    real(dp), dimension(3,3) :: res
    integer(i4) :: i, j
    do j = 1, 3
      do i = 1, 3
        res(i, j) = a(j, i)
      end do
    end do
  end function mat3_transpose

  pure function mat3_vec_mul(a, v) result(res)
    real(dp), dimension(3,3), intent(in) :: a
    real(dp), dimension(3), intent(in) :: v
    real(dp), dimension(3) :: res
    res(1) = a(1,1)*v(1) + a(1,2)*v(2) + a(1,3)*v(3)
    res(2) = a(2,1)*v(1) + a(2,2)*v(2) + a(2,3)*v(3)
    res(3) = a(3,1)*v(1) + a(3,2)*v(2) + a(3,3)*v(3)
  end function mat3_vec_mul

  subroutine mat3_inv(a, a_inv, det, err)
    real(dp), dimension(3,3), intent(in) :: a
    real(dp), dimension(3,3), intent(out) :: a_inv
    real(dp), intent(out) :: det
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: inv_det

    det = mat3_det(a)
    if (abs(det) < 1.0e-30_dp) then
      a_inv = 0.0_dp
      if (present(err)) call err%set(UMBRA_ERR_SINGULAR_METRIC, "Singular 3x3 matrix encountered", "mat3_inv")
      return
    end if

    inv_det = 1.0_dp / det
    a_inv(1,1) =  (a(2,2)*a(3,3) - a(2,3)*a(3,2)) * inv_det
    a_inv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2)) * inv_det
    a_inv(1,3) =  (a(1,2)*a(2,3) - a(1,3)*a(2,2)) * inv_det

    a_inv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1)) * inv_det
    a_inv(2,2) =  (a(1,1)*a(3,3) - a(1,3)*a(3,1)) * inv_det
    a_inv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1)) * inv_det

    a_inv(3,1) =  (a(2,1)*a(3,2) - a(2,2)*a(3,1)) * inv_det
    a_inv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1)) * inv_det
    a_inv(3,3) =  (a(1,1)*a(2,2) - a(1,2)*a(2,1)) * inv_det

    if (present(err)) call err%clear()
  end subroutine mat3_inv

  pure function mat4_det(m) result(det)
    real(dp), dimension(0:3, 0:3), intent(in) :: m
    real(dp) :: det
    real(dp) :: s0, s1, s2, s3, s4, s5
    real(dp) :: c0, c1, c2, c3, c4, c5

    s0 = m(0,0)*m(1,1) - m(0,1)*m(1,0)
    s1 = m(0,0)*m(1,2) - m(0,2)*m(1,0)
    s2 = m(0,0)*m(1,3) - m(0,3)*m(1,0)
    s3 = m(0,1)*m(1,2) - m(0,2)*m(1,1)
    s4 = m(0,1)*m(1,3) - m(0,3)*m(1,1)
    s5 = m(0,2)*m(1,3) - m(0,3)*m(1,2)

    c5 = m(2,2)*m(3,3) - m(2,3)*m(3,2)
    c4 = m(2,1)*m(3,3) - m(2,3)*m(3,1)
    c3 = m(2,1)*m(3,2) - m(2,2)*m(3,1)
    c2 = m(2,0)*m(3,3) - m(2,3)*m(3,0)
    c1 = m(2,0)*m(3,2) - m(2,2)*m(3,0)
    c0 = m(2,0)*m(3,1) - m(2,1)*m(3,0)

    det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0
  end function mat4_det

  pure function mat4_vec_mul(m, v) result(res)
    real(dp), dimension(0:3, 0:3), intent(in) :: m
    real(dp), dimension(0:3), intent(in) :: v
    real(dp), dimension(0:3) :: res
    integer(i4) :: mu, nu
    do mu = 0, 3
      res(mu) = 0.0_dp
      do nu = 0, 3
        res(mu) = res(mu) + m(mu, nu) * v(nu)
      end do
    end do
  end function mat4_vec_mul

  subroutine mat4_inv(m, inv, det, err)
    real(dp), dimension(0:3, 0:3), intent(in) :: m
    real(dp), dimension(0:3, 0:3), intent(out) :: inv
    real(dp), intent(out) :: det
    type(umbra_error_t), intent(inout), optional :: err

    real(dp) :: s0, s1, s2, s3, s4, s5
    real(dp) :: c0, c1, c2, c3, c4, c5
    real(dp) :: inv_det

    s0 = m(0,0)*m(1,1) - m(0,1)*m(1,0)
    s1 = m(0,0)*m(1,2) - m(0,2)*m(1,0)
    s2 = m(0,0)*m(1,3) - m(0,3)*m(1,0)
    s3 = m(0,1)*m(1,2) - m(0,2)*m(1,1)
    s4 = m(0,1)*m(1,3) - m(0,3)*m(1,1)
    s5 = m(0,2)*m(1,3) - m(0,3)*m(1,2)

    c5 = m(2,2)*m(3,3) - m(2,3)*m(3,2)
    c4 = m(2,1)*m(3,3) - m(2,3)*m(3,1)
    c3 = m(2,1)*m(3,2) - m(2,2)*m(3,1)
    c2 = m(2,0)*m(3,3) - m(2,3)*m(3,0)
    c1 = m(2,0)*m(3,2) - m(2,2)*m(3,0)
    c0 = m(2,0)*m(3,1) - m(2,1)*m(3,0)

    det = s0*c5 - s1*c4 + s2*c3 + s3*c2 - s4*c1 + s5*c0
    if (abs(det) < 1.0e-30_dp) then
      inv = 0.0_dp
      if (present(err)) call err%set(UMBRA_ERR_SINGULAR_METRIC, "Singular 4x4 matrix", "mat4_inv")
      return
    end if

    inv_det = 1.0_dp / det

    inv(0,0) = ( m(1,1)*c5 - m(1,2)*c4 + m(1,3)*c3) * inv_det
    inv(0,1) = (-m(0,1)*c5 + m(0,2)*c4 - m(0,3)*c3) * inv_det
    inv(0,2) = ( m(3,1)*s5 - m(3,2)*s4 + m(3,3)*s3) * inv_det
    inv(0,3) = (-m(2,1)*s5 + m(2,2)*s4 - m(2,3)*s3) * inv_det

    inv(1,0) = (-m(1,0)*c5 + m(1,2)*c2 - m(1,3)*c1) * inv_det
    inv(1,1) = ( m(0,0)*c5 - m(0,2)*c2 + m(0,3)*c1) * inv_det
    inv(1,2) = (-m(3,0)*s5 + m(3,2)*s2 - m(3,3)*s1) * inv_det
    inv(1,3) = ( m(2,0)*s5 - m(2,2)*s2 + m(2,3)*s1) * inv_det

    inv(2,0) = ( m(1,0)*c4 - m(1,1)*c2 + m(1,3)*c0) * inv_det
    inv(2,1) = (-m(0,0)*c4 + m(0,1)*c2 - m(0,3)*c0) * inv_det
    inv(2,2) = ( m(3,0)*s4 - m(3,1)*s2 + m(3,3)*s0) * inv_det
    inv(2,3) = (-m(2,0)*s4 + m(2,1)*s2 - m(2,3)*s0) * inv_det

    inv(3,0) = (-m(1,0)*c3 + m(1,1)*c1 - m(1,2)*c0) * inv_det
    inv(3,1) = ( m(0,0)*c3 - m(0,1)*c1 + m(0,2)*c0) * inv_det
    inv(3,2) = (-m(3,0)*s3 + m(3,1)*s1 - m(3,2)*s0) * inv_det
    inv(3,3) = ( m(2,0)*s3 - m(2,1)*s1 + m(2,2)*s0) * inv_det

    if (present(err)) call err%clear()
  end subroutine mat4_inv

  pure subroutine mat3_symm_eigenvalues(a, ev)
    ! Analytic eigenvalues of a real symmetric 3x3 matrix via Cardano trigonometric formula
    real(dp), dimension(3,3), intent(in) :: a
    real(dp), dimension(3), intent(out) :: ev

    real(dp) :: p1, q, p, phi
    real(dp) :: b(3,3)

    p1 = a(1,2)**2 + a(1,3)**2 + a(2,3)**2
    if (p1 < 1.0e-28_dp) then
      ev(1) = a(1,1)
      ev(2) = a(2,2)
      ev(3) = a(3,3)
      return
    end if

    q = (a(1,1) + a(2,2) + a(3,3)) / 3.0_dp
    p1 = (a(1,1) - q)**2 + (a(2,2) - q)**2 + (a(3,3) - q)**2 + 2.0_dp * p1
    p = sqrt(p1 / 6.0_dp)

    b = (a - q * reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3,3])) / p
    phi = mat3_det(b) / 2.0_dp
    if (phi <= -1.0_dp) then
      phi = PI / 3.0_dp
    else if (phi >= 1.0_dp) then
      phi = 0.0_dp
    else
      phi = acos(phi) / 3.0_dp
    end if

    ev(1) = q + 2.0_dp * p * cos(phi)
    ev(3) = q + 2.0_dp * p * cos(phi + TWO_PI / 3.0_dp)
    ev(2) = 3.0_dp * q - ev(1) - ev(3)
  end subroutine mat3_symm_eigenvalues

  subroutine solve_newton_bisection(f_df, x_min, x_max, root, max_iter, tol, err)
    procedure(scalar_func_and_deriv_1d) :: f_df
    real(dp), intent(in) :: x_min, x_max
    real(dp), intent(out) :: root
    integer(i4), intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: iter, max_it
    real(dp) :: target_tol, x, f, df, dx, dx_old
    real(dp) :: xl, xh, fl, fh

    max_it = 80_i4
    if (present(max_iter)) max_it = max_iter
    target_tol = TOL_CONVERGENCE
    if (present(tol)) target_tol = tol

    call f_df(x_min, fl, df)
    call f_df(x_max, fh, df)

    if (fl * fh > 0.0_dp) then
      root = 0.5_dp * (x_min + x_max)
      if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, &
        "Root not bracketed in solve_newton_bisection", "solve_newton_bisection")
      return
    end if

    if (fl < 0.0_dp) then
      xl = x_min
      xh = x_max
    else
      xl = x_max
      xh = x_min
    end if

    root = 0.5_dp * (x_min + x_max)
    dx_old = abs(x_max - x_min)
    dx = dx_old

    call f_df(root, f, df)

    do iter = 1, max_it
      if (((root - xh)*df - f)*((root - xl)*df - f) > 0.0_dp .or. &
          abs(2.0_dp*f) > abs(dx_old*df)) then
        ! Bisect
        dx_old = dx
        dx = 0.5_dp * (xh - xl)
        root = xl + dx
      else
        ! Newton-Raphson step
        dx_old = dx
        dx = f / df
        x = root
        root = root - dx
        if (abs(x - root) <= EPSILON_DP * max(1.0_dp, abs(root))) exit
      end if

      if (abs(dx) < target_tol .or. abs(f) < target_tol) then
        if (present(err)) call err%clear()
        return
      end if

      call f_df(root, f, df)
      if (f < 0.0_dp) then
        xl = root
      else
        xh = root
      end if
    end do

    if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, &
      "Maximum iterations reached in solve_newton_bisection", "solve_newton_bisection")
  end subroutine solve_newton_bisection

  subroutine solve_brent(f, a_in, b_in, root, max_iter, tol, err)
    procedure(scalar_func_1d) :: f
    real(dp), intent(in) :: a_in, b_in
    real(dp), intent(out) :: root
    integer(i4), intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    type(umbra_error_t), intent(inout), optional :: err

    integer(i4) :: iter, max_it
    real(dp) :: target_tol, a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm

    max_it = 100_i4
    if (present(max_iter)) max_it = max_iter
    target_tol = TOL_CONVERGENCE
    if (present(tol)) target_tol = tol

    a = a_in
    b = b_in
    fa = f(a)
    fb = f(b)

    if (fa * fb > 0.0_dp) then
      root = b
      if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, &
        "Root not bracketed in Brent method", "solve_brent")
      return
    end if

    c = b
    fc = fb
    d = b - a
    e = d

    do iter = 1, max_it
      if ((fb > 0.0_dp .and. fc > 0.0_dp) .or. (fb < 0.0_dp .and. fc < 0.0_dp)) then
        c = a
        fc = fa
        d = b - a
        e = d
      end if

      if (abs(fc) < abs(fb)) then
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
      end if

      tol1 = 2.0_dp * EPSILON_DP * abs(b) + 0.5_dp * target_tol
      xm = 0.5_dp * (c - b)

      if (abs(xm) <= tol1 .or. abs(fb) <= target_tol) then
        root = b
        if (present(err)) call err%clear()
        return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb / fa
        if (abs(a - c) < 1.0e-15_dp) then
          p = 2.0_dp * xm * s
          q = 1.0_dp - s
        else
          q = fa / fc
          r = fb / fc
          p = s * (2.0_dp * xm * q * (q - r) - (b - a) * (r - 1.0_dp))
          q = (q - 1.0_dp) * (r - 1.0_dp) * (s - 1.0_dp)
        end if
        if (p > 0.0_dp) q = -q
        p = abs(p)
        if (2.0_dp * p < min(3.0_dp * xm * q - abs(tol1 * q), abs(e * q))) then
          e = d
          d = p / q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if

      a = b
      fa = fb
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1, xm)
      end if
      fb = f(b)
    end do

    root = b
    if (present(err)) call err%set(UMBRA_ERR_CON2PRIM_FAILED, &
      "Maximum iterations reached in Brent method", "solve_brent")
  end subroutine solve_brent

end module umbra_math
