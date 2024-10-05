      subroutine init_phi(ID,ndiv,svall,xl,yl,zl,dx,dy,dz
     & ,bet,rin,xin,yin,zin,phi,nbub,l)

      implicit none
      integer ID,ndiv,svall(3)
      real*8  xl,yl,zl,dx,dy,dz
      real*8 bet,rin,xin,yin,zin,l,nbub
      real*8 phi(-2:svall(1)+3,-2:svall(2)/ndiv+3,-2:svall(3)+3,0:nbub)

      integer ni,nj,nk
      integer i,j,k,jj,n
      real*8 dcell 
      real*8 orgx,orgy,orgz
      real*8 a0,xcent,ycent,zcent
      real*8 x00,y00,z00,xp1,xm1,sdf_bet,q00

      ni=svall(1)
      nj=svall(2)/ndiv
      nk=svall(3)

      dcell=(dx*dy*dz)**(1.0d0/3.0d0)

      orgx=-dx*0.5d0
      orgy=-dy*0.5d0
      orgz=-dz*0.5d0

      a0   =rin                 !initial radius
      xcent=xin
      ycent=yin
      zcent=zin

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,jj)
!$OMP$ PRIVATE(x00,y00,z00,xm1,xp1,sdf_bet,q00)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(dx,dy,dz,orgx,orgy,orgz,xcent,ycent,zcent,xl)
!$OMP$ SHARED(a0,dcell,bet,phi,l)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      jj=j+ID*nj
      x00=dble(i )*dx+orgx-xcent
      y00=dble(jj)*dy+orgy-ycent
      z00=dble(k )*dz+orgz-zcent
      xm1=x00-xl
      xp1=x00+xl
      x00=min(abs(x00),abs(xm1))
      x00=min(    x00 ,abs(xp1))
      sdf_bet=(a0-sqrt(x00**2+y00**2+z00**2))/dcell*bet
      q00=0.0d0
      if(sdf_bet.ge.-5.0d0)then
      if(sdf_bet.le.5.0d0)then
      q00=1.0d0/(1.0d0+exp(-2.0d0*sdf_bet))
      else
      q00=1.0d0
      endif
      phi(i,j,k,l)=max(min(phi(i,j,k,l)+q00,1.0d0),0.0d0)
      endif
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

! lastly phi(i,j,k,0)=sum phi(i,j,k,1:nbub)
      if(l.eq.nbub)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,phi)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      phi(i,j,k,0)=phi(i,j,k,1)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      do n=2,nbub
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(l,ni,nj,nk,nbub,phi,n)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      phi(i,j,k,0)=phi(i,j,k,0)+phi(i,j,k,n)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,nbub,phi)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      phi(i,j,k,0)=min(max(phi(i,j,k,0),0.0d0),1.0d0)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      endif

      return
      end

