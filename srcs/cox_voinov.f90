! cox_voinov.f90
!
! G(θ, λ=1) の 3次近似（30–150°でフィット）
!   c0=0.192131713779, a1=0.230884833059, a3=-0.055730996374
! θ はラジアン。外側角は
!   G(θ_M) = G(θ_m) + S   （S = ±Ca ln(x/L) 等）
! を満たす θ を逆写像で求める。

pure elemental real*8 function cox_voinov(theta_rad) result(G)
  implicit none
  real*8, intent(in) :: theta_rad
  real*8 :: pi, x
  real*8, parameter :: c0 = 0.192131713779d0
  real*8, parameter :: a1 = 0.230884833059d0
  real*8, parameter :: a3 = -0.055730996374d0

  pi = acos(-1.0d0)
  x  = theta_rad - 0.5d0*pi
  G  = c0 + a1*x + a3*x*x*x
end function cox_voinov


pure elemental real*8 function dcox_voinov_dtheta(theta_rad) result(dG)
  implicit none
  real*8, intent(in) :: theta_rad
  real*8 :: pi, x
  real*8, parameter :: a1 = 0.230884833059d0
  real*8, parameter :: a3 = -0.055730996374d0

  pi = acos(-1.0d0)
  x  = theta_rad - 0.5d0*pi
  dG = a1 + 3.0d0*a3*x*x
end function dcox_voinov_dtheta


! Newton で G(θ)=g_target を解く補助関数（30–150°帯での使用を想定）
!  - theta_init: 初期値（ラジアン）
!  - g_target  : 目標の G 値（= G(θ_m) + S）
!  - deg_min/max: 近似の有効帯域（度）
!  - tol, itmax: 収束条件
pure elemental real*8 function cox_voinov_inverse(theta_init, g_target,           &
                                                  deg_min, deg_max, tol, itmax)   &
                                         result(theta_sol)
  implicit none
  real*8, intent(in) :: theta_init, g_target
  real*8, intent(in) :: deg_min, deg_max, tol
  integer, intent(in) :: itmax

  ! 明示インタフェース：呼び出し先が pure elemental であることを示す
  interface
    pure elemental real*8 function cox_voinov(theta_rad) result(G)
      implicit none
      real*8, intent(in) :: theta_rad
    end function cox_voinov
    pure elemental real*8 function dcox_voinov_dtheta(theta_rad) result(dG)
      implicit none
      real*8, intent(in) :: theta_rad
    end function dcox_voinov_dtheta
  end interface

  real*8 :: theta_new, err, pi
  integer :: it

  pi = acos(-1.0d0)
  theta_sol = theta_init

  do it = 1, itmax
    err = cox_voinov(theta_sol) - g_target
    if (abs(err) < tol) exit
    theta_new = theta_sol - err / dcox_voinov_dtheta(theta_sol)

    ! 帯域クランプ（この多項式は 30–150°でフィット）
    if (theta_new < deg_min*pi/180.0d0) theta_new = deg_min*pi/180.0d0
    if (theta_new > deg_max*pi/180.0d0) theta_new = deg_max*pi/180.0d0

    if (abs(theta_new - theta_sol) < tol) then
      theta_sol = theta_new
      exit
    end if
    theta_sol = theta_new
  end do
end function cox_voinov_inverse
