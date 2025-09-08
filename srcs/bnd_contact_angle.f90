!<impose the contact angle boundary condition on qk
!<give the static contact angle theta_0 and grid space dx, dy
subroutine bnd_contact_angle(nID, ni, nj, nk, qk, theta_array, dx, dy, dz)
  implicit none
  include 'param.h'

  integer, intent(in)    :: nID(6), ni, nj, nk
  real*8, intent(inout)  :: qk(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)     :: theta_array(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)     :: dx, dy, dz

  integer :: i, k
  real*8 :: pi, eps, aspect_x, aspect_z
  real*8 :: theta_rad, cos_t, sin_t
  real*8 :: nq_x, nq_y, nq_z, norm

  pi       = acos(-1.0d0)
  eps      = 1.0d-20
  aspect_x = dy / dx
  aspect_z = dy / dz

  !======================================================
  ! Y−面（j = 0 層）: 接触角 BC を x,z 方向にも拡張
  !======================================================
  if (nID(Y_MINUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,cos_t,sin_t, &
!$OMP&                     nq_x,nq_y,nq_z,norm) &
!$OMP&                   SHARED(ni,nj,nk,qk,theta_array,pi,eps,aspect_x,aspect_z)
    do k = -2, nk+3
      do i = -2, ni+3
        !--- 1) 接触角から法線ベクトルを計算 ---
        theta_rad = theta_array(i,1,k) * pi/180.0d0
        cos_t     = cos(theta_rad)
        sin_t     = sin(theta_rad)
        nq_x =  sign(1.0d0, i - ni/2) * sin_t * aspect_x
        nq_z =  sign(1.0d0, i - ni/2) * sin_t * aspect_z
        nq_y =  cos_t
        norm = sqrt(nq_x**2 + nq_y**2 + nq_z**2 + eps)
        nq_x = nq_x / norm
        nq_y = nq_y / norm
        nq_z = nq_z / norm

        !--- 2) 高次外挿: x と z の勾配を考慮 ---
        qk(i,0  ,k) = qk(i,1,k)                                                   &
                     - (nq_y/(nq_x+eps)) * (qk(i+1,1,k) - qk(i-1,1,k)) * 0.5d0*aspect_x &
                     - (nq_y/(nq_z+eps)) * (qk(i,1,k+1) - qk(i,1,k-1)) * 0.5d0*aspect_z

        !--- ゴーストセルを線形反転で設定 ---
        qk(i,-1 ,k) = 2.0d0*qk(i,0  ,k) - qk(i,1,k)
        qk(i,-2 ,k) = 2.0d0*qk(i,-1 ,k) - qk(i,0  ,k)
      end do
    end do
!$OMP END PARALLEL DO
  endif

  !======================================================
  ! Y＋面（j = nj+1,2,3… 層）: 同様に z 方向拡張
  !======================================================
  if (nID(Y_PLUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,cos_t,sin_t, &
!$OMP&                     nq_x,nq_y,nq_z,norm) &
!$OMP&                   SHARED(ni,nj,nk,qk,theta_array,pi,eps,aspect_x,aspect_z)
    do k = -2, nk+3
      do i = -2, ni+3
        theta_rad = theta_array(i,nj,k) * pi/180.0d0
        cos_t     = cos(theta_rad)
        sin_t     = sin(theta_rad)
        nq_x =  sign(1.0d0, i - ni/2) * sin_t * aspect_x
        nq_z =  sign(1.0d0, i - ni/2) * sin_t * aspect_z
        nq_y =  cos_t
        norm = sqrt(nq_x**2 + nq_y**2 + nq_z**2 + eps)
        nq_x = nq_x / norm
        nq_y = nq_y / norm
        nq_z = nq_z / norm

        qk(i,nj+1,k) = qk(i,nj  ,k)                                                   &
                     - (nq_y/(nq_x+eps)) * (qk(i+1,nj,k) - qk(i-1,nj,k)) * 0.5d0*aspect_x &
                     - (nq_y/(nq_z+eps)) * (qk(i,  nj,k+1) - qk(i,  nj,k-1)) * 0.5d0*aspect_z

        qk(i,nj+2,k) = 2.0d0*qk(i,nj+1,k) - qk(i,nj  ,k)
        qk(i,nj+3,k) = 2.0d0*qk(i,nj+2,k) - qk(i,nj+1,k)
      end do
    end do
!$OMP END PARALLEL DO
  endif

  return
end subroutine bnd_contact_angle
