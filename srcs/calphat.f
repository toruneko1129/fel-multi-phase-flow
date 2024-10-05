      subroutine calphat(ni,nj,nk,p,po,phat)
      implicit none
      integer ni,nj,nk
      real*8    p(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   po(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phat(-2:ni+3,-2:nj+3,-2:nk+3)
      
      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(phat,p,po)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-1,ni+3
      phat(i,j,k)=p(i,j,k)*2.0d0-po(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end


