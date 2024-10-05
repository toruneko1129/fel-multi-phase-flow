      subroutine calq(ni,nj,nk,dxinv,dyinv,dzinv,u,v,w,q)
      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8 u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 q(-2:ni+3,-2:nj+3,-2:nk+3)
      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(q,u,v,w,dxinv,dyinv,dzinv)
      do k=-1,nk+3
      do j=-1,nj+3
      do i=-1,ni+3
      q(i,j,k)=
     1  (-u(i-1,j  ,k  )+u(i,j,k))*dxinv
     2 +(-v(i  ,j-1,k  )+v(i,j,k))*dyinv
     2 +(-w(i  ,j  ,k-1)+w(i,j,k))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      return
      end

