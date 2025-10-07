subroutine gnbc(nID, ni, nj, nk, uk, wk, uwall, theta_0_array, &
                surface_tension, zeta_array, theta_array, dx, phi, dbg, delta_x, x_c)

  implicit none
  include 'param.h'
  
  integer nID(6)
  integer ni,nj,nk
  real*8  uk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  wk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  uwall, surface_tension, dx, dy
  real*8  zeta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  theta_0_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  phi(-2:ni+3,-2:nj+3,-2:nk+3)
  
  integer i,k,mi, width
  real*8 theta_0_rad, cos_theta_0, u_cl, cos_theta_t, theta_t_rad, theta_t
  real*8 g_micro, g_macro, coeff, nine_g
  real*8 sgn
  real*8 cos_old
  real*8 u_cl_x, u_cl_z
  real*8 pi, eps, alpha

  integer :: iup_top(-2:nk+3), idn_top(-2:nk+3), iup_bot(-2:nk+3), idn_bot(-2:nk+3)
  real*8  :: xup_top(-2:nk+3), xdn_top(-2:nk+3), xup_bot(-2:nk+3), xdn_bot(-2:nk+3)
  logical :: hup_top(-2:nk+3), hdn_top(-2:nk+3), hup_bot(-2:nk+3), hdn_bot(-2:nk+3)
  logical :: is_band, dbg
  real*8 :: delta_x, x_c

  pi = acos(-1.0d0)
  eps = 1.0d-12
  dy = 13.6d0 / 128

  width = 5  ! 接触点近傍の幅（セル数）
  alpha = 1.0d0
  coeff = 0.917d0 * 1.95d0 * uwall / surface_tension

  call find_two_contacts_on_wall(ni,nj,nk,phi,dx,  1, xup_bot,hup_bot, xdn_bot,hdn_bot)
  call find_two_contacts_on_wall(ni,nj,nk,phi,dx, nj, xup_top,hup_top, xdn_top,hdn_top)
  ! 各壁の2交差を取得（index 版）
  call find_two_contacts_on_wall_index(ni, nj, nk, phi, dx,  1, iup_bot, hup_bot, idn_bot, hdn_bot)
  call find_two_contacts_on_wall_index(ni, nj, nk, phi, dx, nj, iup_top, hup_top, idn_top, hdn_top)
  

  !--- Y−面（j=1）---
  if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,mi,u_cl_x,u_cl_z,u_cl,theta_0_rad,cos_theta_0,cos_theta_t,theta_t_rad,theta_t, is_band, cos_old, g_micro, g_macro, nine_g) &
!$OMP& SHARED(ni,nj,nk,uk,wk,uwall,theta_0_array,surface_tension,zeta_array,theta_array,pi,eps, iup_bot,idn_bot,hup_bot,hdn_bot, width, alpha)
    do k=-2,nk+3
      do i=-2,ni+3
        ! ---- 接触点近傍（±width）チェック：該当しなければスキップ ----
        is_band = .false.
        if (hup_bot(k)) then
          if (abs(i - iup_bot(k)) <= width) is_band = .true.
        end if
        if (hdn_bot(k)) then
          if (abs(i - idn_bot(k)) <= width) is_band = .true.
        end if
        if (.not. is_band) cycle   ! ← 接触点帯域外は更新しない
        ! ------------------------------------------------------------
        mi = mod(i+ni+2, ni+3)

        !--- 接触線速度 x成分 ---
        u_cl_x = (uk(mi,0,k) + uk(mi,1,k) + uk(i,0,k) + uk(i,1,k))/4.0d0 + uwall

        !--- 接触線速度 z成分 ---
        !u_cl_z = (wk(i,0,k) + wk(i,1,k))/2.0d0
        u_cl_z = 0.0d0

        !--- 接触線速度ベクトルの大きさ (符号はx方向で決定) ---
        u_cl = sign(sqrt(u_cl_x**2 + u_cl_z**2 + eps), u_cl_x)
        if (i <= ni/2) then
          u_cl = -u_cl
        endif

        !--- 動的接触角計算 ---
        theta_0_rad = theta_0_array(i, 1, k)*(pi/180.0d0)
        cos_theta_0 = cos(theta_0_rad)
        cos_old = cos(theta_array(i,1,k) * (pi/180.0d0))
        cos_theta_t = cos_theta_0 - zeta_array(i, 1, k)*u_cl/surface_tension
        cos_theta_t = (1.d0 - alpha)*cos_old + alpha*cos_theta_t

        !--- cos_theta_t をクランプ（数値誤差防止） ---
        cos_theta_t = min(1.0d0, max(-1.0d0, cos_theta_t))
        
        theta_t_rad = acos(cos_theta_t)
        !cox-voinov law
        !g_micro = (theta_t_rad**3) / 9.0d0
        !g_macro = g_micro + (0.917d0 * 1.95d0 / surface_tension) * u_cl
        !if (g_macro < 0.d0) then
        !  g_macro = 0.d0   ! 安全側クランプ（数値誤差/極端条件対策）
        !endif
        !nine_g  = 9.0d0 * g_macro
        !if (nine_g <= 0.d0) then
        !  theta_t_rad = 0.d0
        !else
        !  theta_t_rad = nine_g**(1.0d0/3.0d0)
        !end if

        theta_t = theta_t_rad*(180.0d0/pi)

        theta_array(i, 1, k) = theta_t
      enddo
    enddo
!$OMP  END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,k,is_band) SHARED(ni,nj,nk,theta_0_array,theta_array,iup_bot,idn_bot,hup_bot,hdn_bot,width)
do k=-2,nk+3
  do i=-2,ni+3
    is_band = .false.
    if (hup_bot(k)) then
      if (abs(i - iup_bot(k)) <= width) is_band = .true.
    end if
    if (hdn_bot(k)) then
      if (abs(i - idn_bot(k)) <= width) is_band = .true.
    end if
    if (.not. is_band) theta_array(i,1,k) = theta_0_array(i,1,k)
  end do
end do
!$OMP END PARALLEL DO
  endif

  !--- Y＋面（j=nj）---
  if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,mi,u_cl_x,u_cl_z,u_cl,theta_0_rad,cos_theta_0,cos_theta_t,theta_t_rad,theta_t,is_band,cos_old,g_micro,g_macro,nine_g) &
!$OMP& SHARED(ni,nj,nk,uk,wk,uwall,theta_0_array,surface_tension,zeta_array,theta_array,pi,eps, iup_top,idn_top,hup_top,hdn_top, width, alpha,dbg,xup_top,xdn_top,dy,delta_x,x_c)
    do k=-2,nk+3
      do i=-2,ni+3
        ! ---- 接触点近傍（±width）チェック：該当しなければスキップ ----
        is_band = .false.
        if (hup_top(k)) then
          if (abs(i - iup_top(k)) <= width) is_band = .true.
        end if
        if (hdn_top(k)) then
          if (abs(i - idn_top(k)) <= width) is_band = .true.
        end if
        if (.not. is_band) cycle   ! ← 接触点帯域外は更新しない
        ! ------------------------------------------------------------
        mi = mod(i+ni+2, ni+3)

        !--- 接触線速度 x成分 ---
        u_cl_x = (uk(mi,nj,k) + uk(mi,nj+1,k) + uk(i,nj,k) + uk(i,nj+1,k))/4.0d0 - uwall

        !--- 接触線速度 z成分 ---
        !u_cl_z = (wk(i,nj,k) + wk(i,nj+1,k))/2.0d0
        u_cl_z = 0.0d0

        !--- 接触線速度ベクトルの大きさ (符号はx方向で決定) ---
        u_cl = sign(sqrt(u_cl_x**2 + u_cl_z**2 + eps), u_cl_x)
        if (i <= ni/2) then
          u_cl = -u_cl
        endif


        !--- 動的接触角計算 ---
        theta_0_rad = theta_0_array(i, nj, k)*(pi/180.0d0)
        cos_theta_0 = cos(theta_0_rad)
        cos_old = cos(theta_array(i,nj,k) * (pi/180.0d0))
        cos_theta_t = cos_theta_0 - zeta_array(i, nj, k)*u_cl/surface_tension
        cos_theta_t = (1.d0 - alpha)*cos_old + alpha*cos_theta_t

        !--- cos_theta_t をクランプ（数値誤差防止） ---
        cos_theta_t = min(1.0d0, max(-1.0d0, cos_theta_t))
        
        theta_t_rad = acos(cos_theta_t)
        
        !cox-voinov law
        !g_micro = (theta_t_rad**3) / 9.0d0
        !g_macro = g_micro + (0.917d0 * 1.95d0 / surface_tension) * u_cl
        !if (g_macro < 0.d0) then
        !  g_macro = 0.d0   ! 安全側クランプ（数値誤差/極端条件対策）
        !endif
        !nine_g  = 9.0d0 * g_macro
        !if (nine_g <= 0.d0) then
        !  theta_t_rad = 0.d0
        !else
        !  theta_t_rad = nine_g**(1.0d0/3.0d0)
        !end if

        theta_t = theta_t_rad*(180.0d0/pi)

        if (dbg .and. (k==nk/2) .and. (i==idn_top(k))) then
          !write(*,*)  'u=', uwall+u_cl, 'x_c^n=', xdn_top(k) + dy * 2 / tan(theta_array(i, nj, k)*(pi/180.0d0)) - 51.0d0
          delta_x = (uwall+u_cl)*0.25d-2
          x_c = xdn_top(k) + dy * 2 / tan(theta_array(i, nj, k)*(pi/180.0d0)) - 51.0d0
        endif

        theta_array(i, nj, k) = theta_t
      enddo
    enddo
!$OMP  END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,k,is_band) SHARED(ni,nj,nk,theta_0_array,theta_array,iup_top,idn_top,hup_top,hdn_top,width)
do k=-2,nk+3
  do i=-2,ni+3
    is_band = .false.
    if (hup_top(k)) then
      if (abs(i - iup_top(k)) <= width) is_band = .true.
    end if
    if (hdn_top(k)) then
      if (abs(i - idn_top(k)) <= width) is_band = .true.
    end if
    if (.not. is_band) theta_array(i,nj,k) = theta_0_array(i,nj,k)
  end do
end do
!$OMP END PARALLEL DO
  endif

  return
end subroutine gnbc