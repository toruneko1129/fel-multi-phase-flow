      subroutine solphi_mthinc2(ni,nj,nk,dxinv,dyinv,dzinv,dt
     & ,u,v,w,flphix,flphiy,flphiz,phi,phin)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv,dt
      real*8      u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphix(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphiy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphiz(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phin(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k
      real*8 div00

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,div00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,v,w,dxinv,dyinv,dzinv)
!$OMP$ SHARED(phin,phi,dt)
!$OMP$ SHARED(flphix,flphiy,flphiz)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      div00=
     1 +(-u(i-1,j  ,k  )+u(i,j,k))*dxinv
     2 +(-v(i  ,j-1,k  )+v(i,j,k))*dyinv
     3 +(-w(i  ,j  ,k-1)+w(i,j,k))*dzinv 

      phin(i,j,k)=phi(i,j,k)*(1.0d0+dt*div00)
     1 +flphix(i-1,j  ,k  )-flphix(i,j,k)
     2 +flphiy(i  ,j-1,k  )-flphiy(i,j,k)
     3 +flphiz(i  ,j  ,k-1)-flphiz(i,j,k)

      phin(i,j,k)=min(phin(i,j,k),1.0d0)
      phin(i,j,k)=max(phin(i,j,k),0.0d0)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

