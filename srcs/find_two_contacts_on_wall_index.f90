subroutine find_two_contacts_on_wall_index(ni,nj,nk,phi,dx,jwall, i_up,has_up, i_dn,has_dn)
  implicit none
  integer, intent(in)  :: ni, nj, nk, jwall
  real*8,  intent(in)  :: phi(-2:ni+3,-2:nj+3,-2:nk+3), dx
  integer, intent(out) :: i_up(-2:nk+3), i_dn(-2:nk+3)     ! ← 整数で返す
  logical, intent(out) :: has_up(-2:nk+3), has_dn(-2:nk+3)

  integer :: i, k, i_sel
  real*8  :: phi1, phi2, d1, d2

  do k = -2, nk+3
    i_up(k)  = 0
    i_dn(k)  = 0
    has_up(k)= .false.
    has_dn(k)= .false.

    do i = -2, ni+2    ! i+1 を参照するので ni+2 まで
      phi1 = phi(i  , jwall, k)
      phi2 = phi(i+1, jwall, k)

      ! <0.5 → >=0.5（昇順交差）: 左セル i と右セル i+1 の間
      if (.not. has_up(k)) then
        if (phi1 < 0.5d0 .and. phi2 >= 0.5d0) then
          d1 = abs(phi1 - 0.5d0)
          d2 = abs(phi2 - 0.5d0)
          if (d1 <= d2) then
            i_sel = i
          else
            i_sel = i + 1
          end if
          if (i_sel < -2   ) i_sel = -2
          if (i_sel > ni+2 ) i_sel = ni+2
          i_up(k)   = i_sel
          has_up(k) = .true.
        end if
      end if

      ! >0.5 → <=0.5（降順交差）
      if (.not. has_dn(k)) then
        if (phi1 > 0.5d0 .and. phi2 <= 0.5d0) then
          d1 = abs(phi1 - 0.5d0)
          d2 = abs(phi2 - 0.5d0)
          if (d1 <= d2) then
            i_sel = i
          else
            i_sel = i + 1
          end if
          if (i_sel < -2   ) i_sel = -2
          if (i_sel > ni+2 ) i_sel = ni+2
          i_dn(k)   = i_sel
          has_dn(k) = .true.
        end if
      end if

      if (has_up(k) .and. has_dn(k)) exit
    end do
  end do

  !do k = -2, nk+3
  !  if (has_up(k)) then
  !    write(*,'(A,I6,2X,A,I8)') 'find_two_contacts_on_wall_index: k=', k, ' i_up=', i_dn(k)
  !  else
  !    write(*,'(A,I6,2X,A)')    'find_two_contacts_on_wall_index: k=', k, ' i_up: NONE'
  !  end if
  !end do
end subroutine find_two_contacts_on_wall_index
