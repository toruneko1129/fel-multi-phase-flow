subroutine find_interface_positions(ni, nj, nk, phi, dx, dy, xl)
  implicit none
  ! 引数の定義
  integer ni, nj, nk
  real*8 phi(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 dx, dy, xl

  integer :: i, j
  real*8 :: x, y
  real*8 :: phi1, phi2

  ! 界面位置の検出と出力
  x = 0.0d0
  y = (j - 0.5d0) * dy
  do j = 1, nj
    do i = 1, ni
      phi1 = phi(i, j, 1)
      phi2 = phi(i+1, j, 1)
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        ! 線形補間によるx座標の計算
        x += (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
      end if
    end do
  end do
  x /= nj
  print *, x, y

end subroutine find_interface_positions