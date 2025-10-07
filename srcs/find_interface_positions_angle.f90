subroutine find_interface_positions_angle(ni, nj, nk, phi, dx, dy, dz, xl, theta_array, delta_x, x_c, err_gnbc)
  implicit none
  ! 引数の定義
  integer ni, nj, nk
  real*8 phi(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 dx, dy, dz, xl

  integer :: i, j
  real*8 :: x, y
  real*8 :: phi1, phi2
  real*8  theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 pi
  real*8 :: delta_x, x_c, err_gnbc

  pi = acos(-1.0d0)

  ! 界面位置の検出と出力
  do j = nj, nj
    y = (j - 0.5d0) * dy
    do i = 1, ni
      phi1 = phi(i, j, nk/2)
      phi2 = phi(i+1, j, nk/2)
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        ! 線形補間によるx座標の計算
        x = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
		    !write(*,*)  'Δx=', delta_x, 'x_c^(n+1)-x_c^n=', x + dy * 2 / tan(theta_array(i, nj, nk/2)*(pi/180.0d0)) - x_c
        err_gnbc = x + dy * 2 / tan(theta_array(i, nj, nk/2)*(pi/180.0d0)) - x_c - delta_x
        write(*, *) 'err_gnbc=', err_gnbc
        !print *, x, y
      end if
    end do
  end do

end subroutine find_interface_positions_angle