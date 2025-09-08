subroutine init_array_pt_checker(ni,nj,nk,q_a,q_b,q_c,q_array,period,ratio_a)
  implicit none
  integer :: ni,nj,nk
  real*8 :: q_a, q_b, q_c
  real*8 :: q_array(-2:ni+3,-2:nj+3,-2:nk+3)
  integer :: period, ratio_a
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
  !$OMP& SHARED(ni,nj,nk,q_a,q_b,q_c,q_array, period, ratio_a) SCHEDULE(static,1)
  do k = -2, nk+3
    do i = -2, ni+3
      q_array(i, 1, k) = q_a
      if ( mod(i+ni + k+nk, 2) == 0) then
        q_array(i, nj, k) = q_b
      else
        q_array(i, nj, k) = q_c
      endif
      !if ((0 <= k) .AND. (k <= 1) .AND. (0 < i) .AND. (i <= 8) ) then
      !  write(*,'(20e20.10)') q_array(i,nj,k)
      !end if
    end do
  end do
  !$OMP END PARALLEL DO

  return
end subroutine init_array_pt_checker