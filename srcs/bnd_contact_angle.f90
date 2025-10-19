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
  real*8 :: pi, eps, theta_rad, cos_t, sin_t
  real*8 :: gx, gz, gt, nq_x, nq_y, nq_z, norm
  real*8 :: facx, facz, denomx, denomz

  pi   = acos(-1.0d0)
  eps  = 1.0d-20
  facx = 0.5d0 * dy / dx     ! = dy/(2*dx)
  facz = 0.5d0 * dy / dz     ! = dy/(2*dz)

  !======================================================
  ! Y−面（j = 0 層）
  ! 事前に qk の x/z 周期ゴーストが埋まっている前提で
  ! i=1..ni, k=1..nk の中央差分を用いる
  !======================================================
  if (nID(Y_MINUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,cos_t,sin_t,gx,gz,gt, &
!$OMP&                     nq_x,nq_y,nq_z,norm,denomx,denomz) &
!$OMP&                   SHARED(ni,nj,nk,qk,theta_array,pi,eps,facx,facz)
    do k = 1, nk
      do i = 1, ni
        theta_rad = theta_array(i,1,k) * (pi/180.0d0)
        cos_t     = cos(theta_rad)
        sin_t     = sin(theta_rad)

        ! 壁面での接線勾配（中央差分）
        gx = (qk(i+1,1,k)   - qk(i-1,1,k))   / (2.0d0*dx)
        gz = (qk(i,  1,k+1) - qk(i,  1,k-1)) / (2.0d0*dz)
        gt = sqrt(gx*gx + gz*gz)

        ! 法線ベクトル（下壁は +ŷ）
        if (gt .le. eps) then
          nq_x = 0.0d0
          nq_z = 0.0d0
        else
          nq_x = -sin_t * gx / (gt + eps)
          nq_z = -sin_t * gz / (gt + eps)
        end if
        nq_y = +cos_t

        ! 規格化（数値安定）
        norm = sqrt(nq_x*nq_x + nq_y*nq_y + nq_z*nq_z) + eps
        nq_x = nq_x / norm
        nq_y = nq_y / norm
        nq_z = nq_z / norm

        ! 高次外挿：q0 = q1 - (∂q/∂y)*(dy/2)、かつ ∂q/∂y = -(ny/nx)∂x q - (ny/nz)∂z q
        denomx = nq_x
        denomz = nq_z
        qk(i,0,k) = qk(i,1,k)
        if (abs(denomx) .gt. eps) then
          qk(i,0,k) = qk(i,0,k) - (nq_y/(denomx)) * (qk(i+1,1,k)-qk(i-1,1,k)) * facx
        end if
        if (abs(denomz) .gt. eps) then
          qk(i,0,k) = qk(i,0,k) - (nq_y/(denomz)) * (qk(i,1,k+1)-qk(i,1,k-1)) * facz
        end if

        ! ゴーストセル（線形外挿）
        qk(i,-1,k) = 2.0d0*qk(i,0,k)  - qk(i,1,k)
        qk(i,-2,k) = 2.0d0*qk(i,-1,k) - qk(i,0,k)
      end do
    end do
!$OMP END PARALLEL DO
  end if

  !======================================================
  ! Y＋面（j = nj+1, nj+2, …）
  ! 上壁は壁法線が −ŷ である点に注意
  !======================================================
  if (nID(Y_PLUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,cos_t,sin_t,gx,gz,gt, &
!$OMP&                     nq_x,nq_y,nq_z,norm,denomx,denomz) &
!$OMP&                   SHARED(ni,nj,nk,qk,theta_array,pi,eps,facx,facz)
    do k = 1, nk
      do i = 1, ni
        theta_rad = theta_array(i,nj,k) * (pi/180.0d0)
        cos_t     = cos(theta_rad)
        sin_t     = sin(theta_rad)

        gx = (qk(i+1,nj,k)   - qk(i-1,nj,k))   / (2.0d0*dx)
        gz = (qk(i,  nj,k+1) - qk(i,  nj,k-1)) / (2.0d0*dz)
        gt = sqrt(gx*gx + gz*gz)

        ! 法線（上壁は −ŷ）
        if (gt .le. eps) then
          nq_x = 0.0d0
          nq_z = 0.0d0
        else
          nq_x = -sin_t * gx / (gt + eps)
          nq_z = -sin_t * gz / (gt + eps)
        end if
        nq_y = cos_t

        norm = sqrt(nq_x*nq_x + nq_y*nq_y + nq_z*nq_z) + eps
        nq_x = nq_x / norm
        nq_y = nq_y / norm
        nq_z = nq_z / norm

        qk(i,nj+1,k) = qk(i,nj,k)
        denomx = nq_x
        denomz = nq_z
        if (abs(denomx) .gt. eps) then
          qk(i,nj+1,k) = qk(i,nj+1,k) - (nq_y/(denomx)) * (qk(i+1,nj,k)-qk(i-1,nj,k)) * facx
        end if
        if (abs(denomz) .gt. eps) then
          qk(i,nj+1,k) = qk(i,nj+1,k) - (nq_y/(denomz)) * (qk(i,  nj,k+1)-qk(i,  nj,k-1)) * facz
        end if

        qk(i,nj+2,k) = 2.0d0*qk(i,nj+1,k) - qk(i,nj,k)
        qk(i,nj+3,k) = 2.0d0*qk(i,nj+2,k) - qk(i,nj+1,k)
      end do
    end do
!$OMP END PARALLEL DO
  end if

  return
end subroutine bnd_contact_angle
