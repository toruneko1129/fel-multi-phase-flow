subroutine init_theta0(ni,nj,nk,theta_0,theta_0_array)
  implicit none
  integer :: ni,nj,nk
  real*8 :: theta_0
  real*8 :: theta_0_array(-2:ni+3,-2:nj+3,-2:nk+3)
  integer :: i,j,k

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) & 
  !$OMP& SHARED(ni,nj,nk,theta_0_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do j = -2, nj+3
      do i = -2, ni+3
        theta_0_array(i,j,k) = 0.0d0
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,k) & 
  !$OMP& SHARED(ni,nj,nk,theta_0,theta_0_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do i = -2, ni+3
      theta_0_array(i, 1, k) = theta_0
      theta_0_array(i, nj, k) = theta_0
    end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine init_theta0
