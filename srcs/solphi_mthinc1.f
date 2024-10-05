      subroutine solphi_mthinc1(ni,nj,nk
     & ,dxinv,dyinv,dzinv,dt
     & ,bet_mthinc
     & ,phi,phix,phiy,phiz
     & ,u,v,w
     & ,flphix,flphiy,flphiz)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv,dt
      real*8 bet_mthinc
      real*8    phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phix(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phiy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phiz(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphix(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphiy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 flphiz(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k,inr
      real*8 dx,dy,dz
      real*8 bet,gp1,gp2,dcellinv,eps_mthinc,eps_trunc
      real*8 rnx00,rny00,rnz00,q00,phi00
      real*8 coef,xgp1,xgp2,ygp1,ygp2,zgp1,zgp2,wght
      real*8 x1,x2,y1,y2,z1,z2
      real*8 x111,x211,x121,x221,x112,x212,x122,x222
      real*8 ddd,wght_tanh,f,df
      real*8 veldt_dl,xst,xen,yst,yen,zst,zen

      dx=1.0d0/dxinv
      dy=1.0d0/dyinv
      dz=1.0d0/dzinv


!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(flphix,flphiy,flphiz)
      do k=-1,nk+1
      do j=-1,nj+1
      do i=-1,ni+1
      flphix(i,j,k)=0.0d0
      flphiy(i,j,k)=0.0d0
      flphiz(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      bet=bet_mthinc
      gp1=-0.5d0/sqrt(3.0d0)
      gp2= 0.5d0/sqrt(3.0d0)
      dcellinv=1.0d0/sqrt((dx**2+dy**2+dz**2)/3.0d0)
      eps_mthinc=bet/(2.0d0*cosh(-bet*1.5d0)**2)
      eps_trunc=0.5d0*(tanh(-bet*2.5d0)+1.0d0)

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(rnx00,rny00,rnz00,q00,phi00)
!$OMP$ PRIVATE(coef,xgp1,xgp2,ygp1,ygp2,zgp1,zgp2,wght)
!$OMP$ PRIVATE(x1,x2,y1,y2,z1,z2)
!$OMP$ PRIVATE(x111,x211,x121,x221,x112,x212,x122,x222)
!$OMP$ PRIVATE(ddd,inr,wght_tanh,f,df)
!$OMP$ PRIVATE(veldt_dl,xst,xen,yst,yen,zst,zen)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(phix,phiy,phiz,dxinv,dyinv,dzinv)
!$OMP$ SHARED(eps_mthinc,dcellinv,phi,eps_trunc)
!$OMP$ SHARED(u,v,w,flphix,flphiy,flphiz,dt)
!$OMP$ SHARED(gp1,gp2,bet)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1

         rnx00=(
     1 +phix(i-1,j-1,k-1)
     2 +phix(i  ,j-1,k-1)
     3 +phix(i-1,j  ,k-1)
     4 +phix(i  ,j  ,k-1)
     5 +phix(i-1,j-1,k  )
     6 +phix(i  ,j-1,k  )
     7 +phix(i-1,j  ,k  )
     8 +phix(i  ,j  ,k  ))*0.125d0
         rny00=(
     1 +phiy(i-1,j-1,k-1)
     2 +phiy(i  ,j-1,k-1)
     3 +phiy(i-1,j  ,k-1)
     4 +phiy(i  ,j  ,k-1)
     5 +phiy(i-1,j-1,k  )
     6 +phiy(i  ,j-1,k  )
     7 +phiy(i-1,j  ,k  )
     8 +phiy(i  ,j  ,k  ))*0.125d0
         rnz00=(
     1 +phiz(i-1,j-1,k-1)
     2 +phiz(i  ,j-1,k-1)
     3 +phiz(i-1,j  ,k-1)
     4 +phiz(i  ,j  ,k-1)
     5 +phiz(i-1,j-1,k  )
     6 +phiz(i  ,j-1,k  )
     7 +phiz(i-1,j  ,k  )
     8 +phiz(i  ,j  ,k  ))*0.125d0
      q00=sqrt(rnx00**2+rny00**2+rnz00**2)

      if(q00.le.eps_mthinc*dcellinv)then
      phi00=phi(i,j,k)
      if(      phi00.le.eps_trunc)phi00=0.0d0
      if(1.0d0-phi00.le.eps_trunc)phi00=1.0d0

cc
cc<1st-order upwinding
cc
      if(  u(i-1,j,k).le.0.0d0)then
      flphix(i-1,j,k)=u(i-1,j,k)*phi00*dt*dxinv
      endif
      if(  u(i  ,j,k).ge.0.0d0)then
      flphix(i  ,j,k)=u(i  ,j,k)*phi00*dt*dxinv
      endif
      if(  v(i,j-1,k).le.0.0d0)then
      flphiy(i,j-1,k)=v(i,j-1,k)*phi00*dt*dyinv
      endif
      if(  v(i,j  ,k).ge.0.0d0)then
      flphiy(i,j  ,k)=v(i,j  ,k)*phi00*dt*dyinv
      endif
      if(  w(i,j,k-1).le.0.0d0)then
      flphiz(i,j,k-1)=w(i,j,k-1)*phi00*dt*dzinv
      endif
      if(  w(i,j,k  ).ge.0.0d0)then
      flphiz(i,j,k  )=w(i,j,k  )*phi00*dt*dzinv
      endif

      else

cc
cc<mthinc
cc
      coef=1.0d0/q00
      rnx00=rnx00*coef
      rny00=rny00*coef
      rnz00=rnz00*coef
c
c<gauss points
c
      xgp1=gp1
      xgp2=gp2
      ygp1=gp1
      ygp2=gp2
      zgp1=gp1
      zgp2=gp2
      wght=0.125d0
      include'mthinc_x.h'

c
c<find ddd by means of Newton-Raphson method
c
      ddd=0.0d0
      do inr=1,10
      include'mthinc_tanh.h'
      f=wght_tanh-phi(i,j,k)
      df=0.5d0*bet*(
     1 + wght/cosh(bet*(x111+ddd))**2
     2 + wght/cosh(bet*(x211+ddd))**2
     3 + wght/cosh(bet*(x121+ddd))**2
     4 + wght/cosh(bet*(x221+ddd))**2
     5 + wght/cosh(bet*(x112+ddd))**2
     6 + wght/cosh(bet*(x212+ddd))**2
     7 + wght/cosh(bet*(x122+ddd))**2
     8 + wght/cosh(bet*(x222+ddd))**2
     & )
      ddd=ddd-f*df/(df**2+1.0d-99)
      enddo

c
c<compute flux
c

c<xm
      if(u(i-1,j,k).lt.0.0d0)then
      veldt_dl=u(i-1,j,k)*dt*dxinv
      xst=-0.5d0
      xen=-0.5d0-veldt_dl
      xgp1=(xst+xen)*0.5d0+(-xst+xen)*gp1
      xgp2=(xst+xen)*0.5d0+(-xst+xen)*gp2
      ygp1=gp1
      ygp2=gp2
      zgp1=gp1
      zgp2=gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphix(i-1,j,k)=wght_tanh
      endif

c<xp
      if(u(i,j,k).gt.0.0d0)then
      veldt_dl=u(i,j,k)*dt*dxinv
      xst=0.5d0-veldt_dl
      xen=0.5d0
      xgp1=(xst+xen)*0.5d0+(-xst+xen)*gp1
      xgp2=(xst+xen)*0.5d0+(-xst+xen)*gp2
      ygp1=gp1
      ygp2=gp2
      zgp1=gp1
      zgp2=gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphix(i,j,k)=wght_tanh
      endif

c<ym
      if(v(i,j-1,k).lt.0.0d0)then
      veldt_dl=v(i,j-1,k)*dt*dyinv
      xgp1=gp1
      xgp2=gp2
      yst=-0.5d0
      yen=-0.5d0-veldt_dl
      ygp1=(yst+yen)*0.5d0+(-yst+yen)*gp1
      ygp2=(yst+yen)*0.5d0+(-yst+yen)*gp2
      zgp1=gp1
      zgp2=gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphiy(i,j-1,k)=wght_tanh
      endif

c<yp
      if(v(i,j,k).gt.0.0d0)then
      veldt_dl=v(i,j,k)*dt*dyinv
      xgp1=gp1
      xgp2=gp2
      yst=0.5d0-veldt_dl
      yen=0.5d0
      ygp1=(yst+yen)*0.5d0+(-yst+yen)*gp1
      ygp2=(yst+yen)*0.5d0+(-yst+yen)*gp2
      zgp1=gp1
      zgp2=gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphiy(i,j,k)=wght_tanh
      endif

c<zm
      if(w(i,j,k-1).lt.0.0d0)then
      veldt_dl=w(i,j,k-1)*dt*dzinv
      xgp1=gp1
      xgp2=gp2
      ygp1=gp1
      ygp2=gp2
      zst=-0.5d0
      zen=-0.5d0-veldt_dl
      zgp1=(zst+zen)*0.5d0+(-zst+zen)*gp1
      zgp2=(zst+zen)*0.5d0+(-zst+zen)*gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphiz(i,j,k-1)=wght_tanh
      endif

c<zp
      if(w(i,j,k).gt.0.0d0)then
      veldt_dl=w(i,j,k)*dt*dzinv
      xgp1=gp1
      xgp2=gp2
      ygp1=gp1
      ygp2=gp2
      zst=0.5d0-veldt_dl
      zen=0.5d0
      zgp1=(zst+zen)*0.5d0+(-zst+zen)*gp1
      zgp2=(zst+zen)*0.5d0+(-zst+zen)*gp2
      wght=0.125d0*veldt_dl
      include'mthinc_x.h'
      include'mthinc_tanh.h'
      flphiz(i,j,k)=wght_tanh
      endif

      endif

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

