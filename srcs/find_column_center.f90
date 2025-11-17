subroutine find_column_center(ni, nj, nk, phi, dx, dy, dz, xl)
  implicit none
  ! 引数の定義
  integer ni, nj, nk
  real*8 phi(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 dx, dy, dz, xl

  integer :: i, j, k
  real*8 :: x, y, z
  real*8 :: phi1, phi2
  real*8 :: x1_bot, x2_bot, x1_top, x2_top, x_bot, x_top

  ! 界面位置の検出と出力
  do k = 1, nk
    z = (k - 0.5d0) * dz
    do i = 1, ni
      phi1 = phi(i  , nj/2  , k)
      phi2 = phi(i+1, nj/2  , k)
      if (phi1 < 0.5d0 .and. phi2 > 0.5d0) then
        ! 線形補間によるx座標の計算
        x1_bot = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.25d0
      end if
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        ! 線形補間によるx座標の計算
        x2_bot = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
      end if
    end do
    x_bot = x2_bot - x1_bot

    do i = 1, ni
      phi1 = phi(i  , nj/2+1, k)
      phi2 = phi(i+1, nj/2+1, k)
      if (phi1 < 0.5d0 .and. phi2 > 0.5d0) then
        ! 線形補間によるx座標の計算
        x1_top = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.25d0
      end if
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        ! 線形補間によるx座標の計算
        x2_top = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
      end if
    end do
    x_top = x2_top - x1_top
    x = -(x_bot+x_top) / 2
    print *, x, z
  end do
  !end do

end subroutine find_column_center