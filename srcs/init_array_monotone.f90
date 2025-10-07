subroutine init_array_monotone(ni,nj,nk,q_a,q_b,q_c,q_array,period,ratio_a)
  implicit none
  integer :: ni,nj,nk
  real*8 :: q_a, q_b, q_c
  real*8 :: q_array(-2:ni+3,-2:nj+3,-2:nk+3)
  integer :: period, ratio_a
  integer :: i,j,k,ii,kk
  integer :: parity

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
  !$OMP& SHARED(ni,nj,nk,q_a,q_array) SCHEDULE(static,1)
  do k = -2, nk+3
    do i = -2, ni+3
      ! 下側境界
      q_array(i, 1, k) = q_a
      q_array(i, nj, k) = q_a
    end do
  end do
  !$OMP END PARALLEL DO

  ! デバッグ出力（必要なら有効化）
  write(*,'("wall type: monotone")')
  !do k = -2, nk+3
  !  do i = -2, ni+3
  !    if ((4 <= k) .AND. (k <= 5) .AND. (1 <= i) .AND. (i <= 6) ) then
  !      write(*,'(1i2, 1i2, 20e20.10)') i, k, q_array(i,nj,k)
  !    end if
  !  enddo
  !enddo

  return
end subroutine init_array_monotone