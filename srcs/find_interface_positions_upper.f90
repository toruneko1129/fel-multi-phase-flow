subroutine find_interface_positions_upper(ni, nj, nk, phi, dx, dy, dz, xl, theta_array)
  implicit none
  ! 引数の定義
  integer ni, nj, nk
  real*8 phi(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8 dx, dy, dz, xl

  integer :: i, j, k
  real*8 :: x, y, z
  real*8 :: phi1, phi2
  real*8  theta_array(-2:ni+3,-2:nj+3,-2:nk+3)
  integer :: idx

  ! 界面位置の検出と出力
  do k = 1, nk
    z = (k - 0.5d0) * dz
    do i = 1, ni
      phi1 = phi(i, nj, k)
      phi2 = phi(i+1, nj, k)
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        ! 線形補間によるx座標の計算
        x = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
        print *, x, z
      end if
    end do
  end do

  write(*, *) '--- Contact Angle Distribution on Upper Wall ---'

  ! 接触角分布の出力
  do k = 1, nk
    z = (k - 0.5d0) * dz
    do i = 1, ni
      phi1 = phi(i, nj, k)
      phi2 = phi(i+1, nj, k)
      if (phi1 > 0.5d0 .and. phi2 < 0.5d0) then
        x = (i - 0.5d0) * dx + (0.5d0 - phi1) * dx / (phi2 - phi1) - xl * 0.75d0
        idx = ni*3/4 + 1 + INT( x/dx )
        print *, theta_array(idx,nj,k), z
      end if
    end do
  end do

end subroutine find_interface_positions_upper