!=========================================================
!  CSF用 VOF 平滑化
!   - x,z: periodic (bnd_periodic が両方を埋める前提)
!   - y  : 接触角 theta を満たすよう ghost(phi) を再構成 (bnd_contact_angle)
!   - 平滑化は内点(i=1..ni,j=1..nj,k=1..nk) で計算
!=========================================================
subroutine smooth_vof_csf(nID, ni, nj, nk, phi, phihat, theta_array, dx, dy, dz, niter)
  implicit none
  include 'param.h'

  integer, intent(in) :: nID(6), ni, nj, nk, niter
  real*8, intent(in)  :: phi       (-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(out) :: phihat    (-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)  :: theta_array(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)  :: dx, dy, dz

  integer :: it, i, j, k
  real*8  :: tmp(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8  :: w0, wn, val
  logical :: is3d

  is3d = (nk >= 3)

  ! 重み（2D: center 3/4 + 4近傍 1/16, 3D: center 3/4 + 6近傍 1/24）
  w0 = 0.75d0
  if (is3d) then
    wn = (1.0d0 - w0) / 6.0d0
  else
    wn = (1.0d0 - w0) / 4.0d0
  end if

  call cpy(ni, nj, nk, phi, tmp)

  ! ★最初に tmp のゴーストも整合させる（順序：周期→接触角→周期）
  call bnd_periodic(ni, nj, nk, tmp)
  call bnd_contact_angle(nID, ni, nj, nk, tmp, theta_array, dx, dy, dz)
  !call bnd_neumann(nID,ni,nj,nk,tmp)
  call bnd_periodic(ni, nj, nk, tmp)

  do it = 1, niter
    call cpy(ni, nj, nk, tmp, phihat)

    if (is3d) then
!$OMP PARALLEL DO PRIVATE(i,j,val) &
!$OMP& SHARED(ni,nj,nk,tmp,phihat,w0,wn)
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            val = w0*tmp(i,j,k) + wn*( &
                  tmp(i+1,j,k) + tmp(i-1,j,k) + &
                  tmp(i,j+1,k) + tmp(i,j-1,k) + &
                  tmp(i,j,k+1) + tmp(i,j,k-1) )
            if (val < 0.0d0) val = 0.0d0
            if (val > 1.0d0) val = 1.0d0
            phihat(i,j,k) = val
          end do
        end do
      end do
!$OMP END PARALLEL DO
    else
      ! nk=1 or 2 を「2Dフィルタ」で扱う（k方向は触らない）
!$OMP PARALLEL DO PRIVATE(i,j,val) &
!$OMP& SHARED(ni,nj,nk,tmp,phihat,w0,wn)
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            val = w0*tmp(i,j,k) + wn*( &
                  tmp(i+1,j,k) + tmp(i-1,j,k) + &
                  tmp(i,j+1,k) + tmp(i,j-1,k) )
            if (val < 0.0d0) val = 0.0d0
            if (val > 1.0d0) val = 1.0d0
            phihat(i,j,k) = val
          end do
        end do
      end do
!$OMP END PARALLEL DO
    end if

    ! ★毎反復：境界を必ず再構成（順序：周期→接触角→周期）
    call bnd_periodic(ni, nj, nk, phihat)
    call bnd_contact_angle(nID, ni, nj, nk, phihat, theta_array, dx, dy, dz)
	!call bnd_neumann(nID,ni,nj,nk,phihat)
    call bnd_periodic(ni, nj, nk, phihat)

    call cpy(ni, nj, nk, phihat, tmp)
  end do

end subroutine smooth_vof_csf