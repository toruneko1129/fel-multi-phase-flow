!<impose the contact angle boundary condition on qk
!<give the static contact angle theta_0 and grid space dx, dy
subroutine bnd_contact_angle(nID, ni, nj, nk, qk, theta_array, dx, dy)
    
  implicit none
  include 'param.h'
  integer nID(6), ni, nj, nk
  real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 dx, dy

  integer i, k
  real*8 theta_rad, cos_theta, sin_theta, eps, aspect
  real*8 nq_x, nq_y, nq_z, norm

  eps = 1.0d-20
  aspect = dy / dx

  if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,theta_rad,sin_theta,cos_theta,nq_x,nq_y,nq_z,norm) &
!$OMP& SHARED(ni,nj,nk,qk,theta_array,eps,aspect)
    do k=-2,nk+3
      do i=-2,ni+3
        theta_rad = theta_array(i,1,k) * (acos(-1.0d0) / 180.0d0)
        cos_theta = cos(theta_rad)
        sin_theta = sin(theta_rad)

        if (i <= ni / 2) then
          nq_x = -sin_theta * aspect
        else
          nq_x = sin_theta * aspect
        end if
        nq_y = cos_theta
        norm = sqrt(nq_x**2 + nq_y**2 + 1.0d-20)
        nq_x = nq_x / norm
        nq_y = nq_y / norm

        ! 壁面でのqkの値を高次の外挿法で計算
        qk(i,0,k) = qk(i,1,k) - (nq_y / (nq_x + eps)) * (qk(i+1,1,k) - qk(i-1,1,k)) * 0.5d0 * aspect
        ! ゴーストセルの値を設定（必要に応じて）
        qk(i,-1,k) = 2.0d0 * qk(i,0,k) - qk(i,1,k)
        qk(i,-2,k) = 2.0d0 * qk(i,-1,k) - qk(i,0,k)
      enddo
    enddo
!$OMP  END PARALLEL DO
  endif

  if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO PRIVATE(i,k,theta_rad,sin_theta,cos_theta,nq_x,nq_y,nq_z,norm) &
!$OMP& SHARED(ni,nj,nk,qk,theta_array,eps,aspect)
    do k=-2,nk+3
      do i=-2,ni+3
        theta_rad = theta_array(i,nj,k) * (acos(-1.0d0) / 180.0d0)
        cos_theta = cos(theta_rad)
        sin_theta = sin(theta_rad)

        if (i < ni / 2) then
          nq_x = -sin_theta * aspect
        else
          nq_x = sin_theta * aspect
        end if
        nq_y = cos_theta
        norm = sqrt(nq_x**2 + nq_y**2 + 1.0d-20)
        nq_x = nq_x / norm
        nq_y = nq_y / norm

        ! 壁面でのqkの値を高次の外挿法で計算
        qk(i,nj+1,k) = qk(i,nj,k) - (nq_y / (nq_x + eps)) * (qk(i+1,nj,k) - qk(i-1,nj,k)) * 0.5d0 * aspect
        ! ゴーストセルの値を設定（必要に応じて）
        qk(i,nj+2,k) = 2.0d0 * qk(i,nj+1,k) - qk(i,nj,k)
        qk(i,nj+3,k) = 2.0d0 * qk(i,nj+2,k) - qk(i,nj+1,k)
      enddo
    enddo
!$OMP  END PARALLEL DO
  endif

  return
end
