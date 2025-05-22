subroutine init_array(ni,nj,nk,q,q_array)
  implicit none
  integer :: ni,nj,nk
  real*8 :: q
  real*8 :: q_array(-2:ni+3,-2:nj+3,-2:nk+3)
  integer :: i,j,k

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) & 
  !$OMP& SHARED(ni,nj,nk,q_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do j = -2, nj+3
      do i = -2, ni+3
        q_array(i,j,k) = 0.0d0
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,k) & 
  !$OMP& SHARED(ni,nj,nk,q,q_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do i = -2, ni+3
      q_array(i, 1, k) = q
      q_array(i, nj, k) = q
    end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine init_array
