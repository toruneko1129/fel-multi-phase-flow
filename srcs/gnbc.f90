subroutine gnbc(nID, ni, nj, nk, uk, uwall, theta_0, surface_tension, zeta, theta_array)

  implicit none
  include 'param.h'
  integer nID(6)
  integer ni,nj,nk
  real*8  uk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  uwall, theta_0, surface_tension, zeta
  real*8  theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  theta_0_rad, cos_theta_0, u_cl, cos_theta_t, theta_t_rad, theta_t

  integer i,k

  theta_0_rad = theta_0 * (acos(-1.0d0) / 180.0d0)
  cos_theta_0 = cos(theta_0_rad)

  if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,uk,u_cl,cos_theta_t,theta_t_rad,theta_t) &
!$OMP& SHARED(ni,nj,nk,uwall,theta_0,surface_tension,zeta,theta_array,cos_theta_0)
    do k=-2,nk+3
      do i=-2,ni+3
        u_cl = (uk(mod(i+ni+2,ni+3), 0, k) + uk(mod(i+ni+2,ni+3), 1, k) + &
                 uk(i, 0, k) + uk(i, 1, k)) / 4.0d0 + uwall
        if (i <= ni/2) then
          u_cl = -u_cl
        endif
        cos_theta_t = cos_theta_0 - zeta * u_cl / surface_tension
        theta_t_rad = acos(cos_theta_t)
        theta_t = theta_t_rad * (180.0d0 / acos(-1.0d0))

        theta_array(i, 1, k) = theta_t
      enddo
    enddo
  endif

  if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,uk,u_cl,cos_theta_t,theta_t_rad,theta_t) &
!$OMP& SHARED(ni,nj,nk,uwall,theta_0,surface_tension,zeta,theta_array,cos_theta_0)
    do k=-2,nk+3
      do i=-2,ni+3
        u_cl = (uk(mod(i+ni+2,ni+3), nj, k) + uk(mod(i+ni+2,ni+3), nj+1, k) + &
                 uk(i, 0, k) + uk(i, 1, k)) / 4.0d0 - uwall
        if (i <= ni/2) then
          u_cl = -u_cl
        endif
        cos_theta_t = cos_theta_0 - zeta * u_cl / surface_tension
        theta_t_rad = acos(cos_theta_t)
        theta_t = theta_t_rad * (180.0d0 / acos(-1.0d0))

        theta_array(i, nj, k) = theta_t
      enddo
    enddo
  endif

  return
end subroutine gnbc