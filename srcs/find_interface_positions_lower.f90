subroutine find_interface_positions_lower(ni, nj, nk, phi, dx, dy, dz, xl)
  implicit none

  integer, intent(in)    :: ni, nj, nk
  real*8, intent(in)     :: phi(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)     :: dx, dy, dz, xl

  integer :: i, j, k
  real*8 :: sum_x, x_k, x_avg, y
  real*8 :: phi1, phi2

  sum_x = 0.0d0

  ! 各 k 平面で界面位置 x_k を計算
  do k = 1, nk
    x_k = 0.0d0
    ! しきい値を超える最初の i で線形補間
    do i = 1, ni
      phi1 = phi(i,   1, k)
      phi2 = phi(i+1, 1, k)
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        x_k = (i - 0.5d0)*dx + (0.5d0 - phi1)*dx/(phi2 - phi1) - xl*0.75d0
        exit  ! 各 k につき１回検出できればループを抜ける
      end if
    end do
    sum_x = sum_x + x_k
  end do

  ! k 平面での平均を取って出力
  x_avg = sum_x / nk
  print *, 'Lower Interface Position (Average over k):'
  print *, -x_avg


end subroutine find_interface_positions_lower