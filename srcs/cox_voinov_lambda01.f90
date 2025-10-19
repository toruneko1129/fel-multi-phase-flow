! cox_voinov_lambda01.f90
!
! G(θ, λ=0.1) の 4次近似（30–150°でフィット）
!   c0=0.355, a1=0.583, a2=0.201, a3=-0.120, a4=-0.0728
! θ はラジアン。外側角は
!   G(θ_M) = G(θ_m) + S   （S = ±Ca ln(x/L) 等）
! を満たす θ を逆写像で求める。

pure elemental real*8 function cox_voinov_lambda01(theta_rad) result(G)
  implicit none
  real*8, intent(in) :: theta_rad
  real*8 :: pi, x
  real*8, parameter :: c0 = 0.355d0
  real*8, parameter :: a1 = 0.583d0
  real*8, parameter :: a2 = 0.201d0
  real*8, parameter :: a3 = -0.120d0
  real*8, parameter :: a4 = -0.0728d0

  pi = acos(-1.0d0)
  x  = theta_rad - 0.5d0*pi
  G  = ((a4*x + a3)*x + a2)*x*x + a1*x + c0
  ! = a4 x^4 + a3 x^3 + a2 x^2 + a1 x + c0
end function cox_voinov_lambda01


pure elemental real*8 function dcox_voinov_dtheta_lambda01(theta_rad) result(dG)
  implicit none
  real*8, intent(in) :: theta_rad
  real*8 :: pi, x
  real*8, parameter :: a1 = 0.583d0
  real*8, parameter :: a2 = 0.201d0
  real*8, parameter :: a3 = -0.120d0
  real*8, parameter :: a4 = -0.0728d0

  pi = acos(-1.0d0)
  x  = theta_rad - 0.5d0*pi
  dG = a1 + 2.0d0*a2*x + 3.0d0*a3*x*x + 4.0d0*a4*x*x*x
end function dcox_voinov_dtheta_lambda01


! Newton で G(θ)=g_target を解く補助関数（30–150°帯での使用を想定）
!  - theta_init: 初期値（ラジアン）
!  - g_target  : 目標の G 値（= G(θ_m) + S）
!  - deg_min/max: 近似の有効帯域（度）
!  - tol, itmax: 収束条件
pure elemental real*8 function cox_voinov_inverse_lambda01(theta_init, g_target,           &
                                                   deg_min, deg_max, tol, itmax)  &
                                          result(theta_sol)
  implicit none
  real*8, intent(in) :: theta_init, g_target
  real*8, intent(in) :: deg_min, deg_max, tol
  integer, intent(in) :: itmax

  interface
    pure elemental real*8 function cox_voinov_lambda01(theta_rad) result(G)
      implicit none
      real*8, intent(in) :: theta_rad
    end function cox_voinov_lambda01
    pure elemental real*8 function dcox_voinov_dtheta_lambda01(theta_rad) result(dG)
      implicit none
      real*8, intent(in) :: theta_rad
    end function dcox_voinov_dtheta_lambda01
  end interface

  real*8 :: theta_new, err, pi, th_min, th_max
  integer :: it

  pi = acos(-1.0d0)
  th_min = deg_min*pi/180.0d0
  th_max = deg_max*pi/180.0d0

  theta_sol = theta_init

  do it = 1, itmax
    err = cox_voinov_lambda01(theta_sol) - g_target
    if (abs(err) < tol) exit

    theta_new = theta_sol - err / dcox_voinov_dtheta_lambda01(theta_sol)

    ! 帯域クランプ（この多項式は 30–150°でフィット）
    if (theta_new < th_min) theta_new = th_min
    if (theta_new > th_max) theta_new = th_max

    if (abs(theta_new - theta_sol) < tol) then
      theta_sol = theta_new
      exit
    end if
    theta_sol = theta_new
  end do
end function cox_voinov_inverse_lambda01
