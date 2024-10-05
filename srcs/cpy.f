      subroutine cpy(ni,nj,nk,q1,q2)
      implicit none
      integer ni,nj,nk
      real*8 q1(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 q2(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(q2,q1)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      q2(i,j,k)=q1(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
