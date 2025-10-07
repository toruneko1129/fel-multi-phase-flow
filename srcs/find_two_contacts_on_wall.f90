subroutine find_two_contacts_on_wall(ni,nj,nk,phi,dx,jwall, x_up,has_up, x_dn,has_dn)
  implicit none
  integer, intent(in) :: ni, nj, nk, jwall
  real*8,  intent(in) :: phi(-2:ni+3,-2:nj+3,-2:nk+3), dx
  real*8,  intent(out):: x_up(-2:nk+3), x_dn(-2:nk+3)
  logical, intent(out):: has_up(-2:nk+3), has_dn(-2:nk+3)
  integer :: i, k
  real*8 :: phi1, phi2, x0, eps
  eps = 1.0d-14  ! 0割回避

  do k=-2, nk+3
    x_up(k)=0.0d0; x_dn(k)=0.0d0
    has_up(k)=.false.; has_dn(k)=.false.
    do i=-2, ni+2                      ! i+1 を読むので ni+2 まで
      phi1 = phi(i  , jwall, k)
      phi2 = phi(i+1, jwall, k)

      ! <0.5 → >=0.5（昇順交差）
      if (.not.has_up(k)) then
        if ( (phi1<0.5d0 .and. phi2>=0.5d0)) then
          x0 = (dble(i)-0.5d0)*dx + (0.5d0-phi1)*dx/(phi2-phi1)
          x_up(k)=x0; has_up(k)=.true.
        end if
      end if

      ! >0.5 → <=0.5（降順交差）
      if (.not.has_dn(k)) then
        if ( (phi1>0.5d0 .and. phi2<=0.5d0)) then
          x0 = (dble(i)-0.5d0)*dx + (0.5d0-phi1)*dx/(phi2-phi1)
          x_dn(k)=x0; has_dn(k)=.true.
        end if
      end if

      if (has_up(k) .and. has_dn(k)) exit
    end do
  end do

  ! do k = -2, nk+3
  !   if (has_up(k)) then
  !     write(*,'(A,I6,2X,A,I8)') 'find_two_contacts_on_wall_index: k=', k, ' i_up=', i_up(k)
  !   else
  !     write(*,'(A,I6,2X,A)')    'find_two_contacts_on_wall_index: k=', k, ' i_up: NONE'
  !   end if
  ! end do
end subroutine find_two_contacts_on_wall
