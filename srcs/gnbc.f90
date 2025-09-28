subroutine gnbc(nID, ni, nj, nk, uk, wk, uwall, theta_0_array, &
                surface_tension, zeta_array, theta_array)

  implicit none
  include 'param.h'
  
  integer nID(6)
  integer ni,nj,nk
  real*8  uk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  wk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  uwall, surface_tension
  real*8  zeta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  theta_0_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  
  integer i,k,mi
  real*8 theta_0_rad, cos_theta_0, u_cl, cos_theta_t, theta_t_rad, theta_t
  real*8 u_cl_x, u_cl_z
  real*8 pi, eps

  pi = acos(-1.0d0)
  eps = 1.0d-12

  !--- Y−面（j=0）---
  if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,mi,u_cl_x,u_cl_z,u_cl,theta_0_rad,cos_theta_0,cos_theta_t,theta_t_rad,theta_t) &
!$OMP& SHARED(ni,nj,nk,uk,wk,uwall,theta_0_array,surface_tension,zeta_array,theta_array,pi,eps)
    do k=-2,nk+3
      do i=-2,ni+3
        mi = mod(i+ni+2, ni+3)

        !--- 接触線速度 x成分 ---
        u_cl_x = (uk(mi,0,k) + uk(mi,1,k) + uk(i,0,k) + uk(i,1,k))/4.0d0 + uwall

        !--- 接触線速度 z成分 ---
        u_cl_z = (wk(mi,0,k) + wk(mi,1,k) + wk(i,0,k) + wk(i,1,k))/4.0d0

        !--- 接触線速度ベクトルの大きさ (符号はx方向で決定) ---
        u_cl = sign(sqrt(u_cl_x**2 + u_cl_z**2 + eps), u_cl_x)
        if (i <= ni/2) then
          u_cl = -u_cl
        endif

        !--- 動的接触角計算 ---
        theta_0_rad = theta_0_array(i, 1, k)*(pi/180.0d0)
        cos_theta_0 = cos(theta_0_rad)
        cos_theta_t = cos_theta_0 - zeta_array(i, 1, k)*u_cl/surface_tension

        !--- cos_theta_t をクランプ（数値誤差防止） ---
        cos_theta_t = min(1.0d0, max(-1.0d0, cos_theta_t))
        
        theta_t_rad = acos(cos_theta_t)
        theta_t = theta_t_rad*(180.0d0/pi)

        theta_array(i, 1, k) = theta_t
      enddo
    enddo
!$OMP  END PARALLEL DO
  endif

  !--- Y＋面（j=nj）---
  if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,mi,u_cl_x,u_cl_z,u_cl,theta_0_rad,cos_theta_0,cos_theta_t,theta_t_rad,theta_t) &
!$OMP& SHARED(ni,nj,nk,uk,wk,uwall,theta_0_array,surface_tension,zeta_array,theta_array,pi,eps)
    do k=-2,nk+3
      do i=-2,ni+3
        mi = mod(i+ni+2, ni+3)

        !--- 接触線速度 x成分 ---
        u_cl_x = (uk(mi,nj,k) + uk(mi,nj+1,k) + uk(i,nj,k) + uk(i,nj+1,k))/4.0d0 - uwall

        !--- 接触線速度 z成分 ---
        u_cl_z = (wk(mi,nj,k) + wk(mi,nj+1,k) + wk(i,nj,k) + wk(i,nj+1,k))/4.0d0

        !--- 接触線速度ベクトルの大きさ (符号はx方向で決定) ---
        u_cl = sign(sqrt(u_cl_x**2 + u_cl_z**2 + eps), u_cl_x)
        if (i <= ni/2) then
          u_cl = -u_cl
        endif


        !--- 動的接触角計算 ---
        theta_0_rad = theta_0_array(i, nj, k)*(pi/180.0d0)
        cos_theta_0 = cos(theta_0_rad)
        cos_theta_t = cos_theta_0 - zeta_array(i, nj, k)*u_cl/surface_tension

        !--- cos_theta_t をクランプ（数値誤差防止） ---
        cos_theta_t = min(1.0d0, max(-1.0d0, cos_theta_t))
        
        theta_t_rad = acos(cos_theta_t)
        theta_t = theta_t_rad*(180.0d0/pi)

        theta_array(i, nj, k) = theta_t
      enddo
    enddo
!$OMP  END PARALLEL DO
  endif

  return
end subroutine gnbc