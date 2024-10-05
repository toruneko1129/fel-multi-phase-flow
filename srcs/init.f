      subroutine init(ni,nj,nk,u,v,w,p,uo,vo,wo,po,dp,phi,nbub)

      implicit none
      integer ni,nj,nk,nbub
      real*8   u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   p(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  uo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  vo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  wo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  po(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  dp(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phi(-2:ni+3,-2:nj+3,-2:nk+3,0:nbub)

      integer i,j,k,l

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,l)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,v,w,p,uo,vo,wo,po,dp,phi,nbub)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
        u(i,j,k)=0.0d0
        v(i,j,k)=0.0d0
        w(i,j,k)=0.0d0
        p(i,j,k)=0.0d0
       uo(i,j,k)=0.0d0
       vo(i,j,k)=0.0d0
       wo(i,j,k)=0.0d0
       po(i,j,k)=0.0d0
       dp(i,j,k)=0.0d0
      do l=0,nbub
       phi(i,j,k,l)=0.0d0
      enddo
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
