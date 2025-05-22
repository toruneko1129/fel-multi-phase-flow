subroutine init_array(ni,nj,nk,q_a,q_b,q_array)
  implicit none
  integer :: ni,nj,nk
  real*8 :: q_a, q_b
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
  !$OMP& SHARED(ni,nj,nk,q_a,q_b,q_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do i = -2, ni+3
      if ( mod(i,2) == 1 ) then
        q_array(i, 1, k) = q_a
        q_array(i, nj, k) = q_a
      else
        q_array(i, 1, k) = q_b
        q_array(i, nj, k) = q_b
      endif
    end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine init_array
