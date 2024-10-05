      subroutine corunp_explicit(nID,ni,nj,nk,rho,dxinv,dyinv,dzinv,dt
     & ,dp,phat,un,vn,wn,pn)

      implicit none
      include 'param.h'
      integer nID(6)
      integer ni,nj,nk
      real*8 rho,dxinv,dyinv,dzinv,dt
      real*8   dp(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phat(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   un(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   wn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   pn(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k,jen

      jen=nj
      if(nID(Y_PLUS).lt.0)jen=nj-1

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(un,dp,rho,dxinv,dt)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      un(i,j,k)=un(i,j,k)+(dp(i,j,k)-dp(i+1,j,k))/rho*dxinv*dt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,jen,nk)
!$OMP$ SHARED(vn,dp,rho,dyinv,dt)
      do k=1,nk
      do j=1,jen
      do i=1,ni
      vn(i,j,k)=vn(i,j,k)+(dp(i,j,k)-dp(i,j+1,k))/rho*dyinv*dt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(wn,dp,rho,dzinv,dt)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      wn(i,j,k)=wn(i,j,k)+(dp(i,j,k)-dp(i,j,k+1))/rho*dzinv*dt
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(pn,phat,dp)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      pn(i,j,k)=phat(i,j,k)+dp(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
