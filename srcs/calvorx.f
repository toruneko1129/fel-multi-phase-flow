      subroutine calvorx(ni,nj,nk,dyinv,dzinv,v,w,vorx)
      implicit none
      integer ni,nj,nk
      real*8 dyinv,dzinv
      real*8    v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 vorx(-2:ni+3,-2:nj+3,-2:nk+3)
      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(vorx,w,v,dyinv,dzinv)
      do k=-2,nk+2
      do j=-2,nj+2
      do i=-2,ni+3
       vorx(i,j,k)=
     1  (-w(i,j,k)+w(i,j+1,k  ))*dyinv
     2 -(-v(i,j,k)+v(i,j  ,k+1))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
