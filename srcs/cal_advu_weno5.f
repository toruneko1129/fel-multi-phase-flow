      subroutine cal_advu_weno5(ni,nj,nk,dxinv,dyinv,dzinv
     & ,u,v,w,adv_u,adv_v,adv_w)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8      u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

      real*8 verysmall,d1_12
      real*8 ucm1,ucp1,vcm1,vcp1,wcm1,wcp1
      real*8 dfm5,dfm3,dfm1,dfp1,dfp3,dfp5
      real*8 vel00
      real*8 dfm_1,dfm_2,dfm_3
      real*8 dfp_1,dfp_2,dfp_3
      real*8 sm_1,sm_2,sm_3
      real*8 sp_1,sp_2,sp_3
      real*8 am_1,am_2,am_3
      real*8 ap_1,ap_2,ap_3
      real*8 deninv_am,deninv_ap
      real*8 wgm_1,wgm_3,wgm_2
      real*8 wgp_1,wgp_2,wgp_3
      real*8 velab,velm1,velp1,adv_weno5

      verysmall=1.0d-80
      d1_12=1.0d0/12.0d0

ccc
ccc<initially clear the terms
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(adv_u,adv_v,adv_w)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      adv_u(i,j,k)=0.0d0
      adv_v(i,j,k)=0.0d0
      adv_w(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
 

c
c<u
c
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(dfm5,dfm3,dfm1,dfp1,dfp3,dfp5)
!$OMP$ PRIVATE(vel00)
!$OMP$ PRIVATE(dfm_1,dfm_2,dfm_3)
!$OMP$ PRIVATE(dfp_1,dfp_2,dfp_3)
!$OMP$ PRIVATE(sm_1,sm_2,sm_3)
!$OMP$ PRIVATE(sp_1,sp_2,sp_3)
!$OMP$ PRIVATE(am_1,am_2,am_3)
!$OMP$ PRIVATE(ap_1,ap_2,ap_3)
!$OMP$ PRIVATE(deninv_am,deninv_ap)
!$OMP$ PRIVATE(wgm_1,wgm_3,wgm_2)
!$OMP$ PRIVATE(wgp_1,wgp_2,wgp_3)
!$OMP$ PRIVATE(velab,velm1,velp1,adv_weno5)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
!$OMP$ SHARED(u,v,w,adv_u)
!$OMP$ SHARED(verysmall,d1_12)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      ucm1=(u(i-1,j  ,k  )+u(i  ,j  ,k  ))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i+1,j  ,k  ))*0.5d0
      vel00=(ucm1+ucp1)*0.5d0
      dfm5= -u(i-3,j,k)+u(i-2,j,k)
      dfm3= -u(i-2,j,k)+u(i-1,j,k)
      dfm1= -u(i-1,j,k)+u(i  ,j,k)
      dfp1= -u(i  ,j,k)+u(i+1,j,k)
      dfp3= -u(i+1,j,k)+u(i+2,j,k)
      dfp5= -u(i+2,j,k)+u(i+3,j,k)
      include'weno5.h'
      adv_u(i,j,k)=adv_weno5*dxinv
 
      vcm1=(v(i  ,j-1,k  )+v(i+1,j-1,k  ))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i+1,j  ,k  ))*0.5d0
      vel00=(vcm1+vcp1)*0.5d0
      dfm5= -u(i,j-3,k)+u(i,j-2,k)
      dfm3= -u(i,j-2,k)+u(i,j-1,k)
      dfm1= -u(i,j-1,k)+u(i,j  ,k)
      dfp1= -u(i,j  ,k)+u(i,j+1,k)
      dfp3= -u(i,j+1,k)+u(i,j+2,k)
      dfp5= -u(i,j+2,k)+u(i,j+3,k)
      include'weno5.h'
      adv_u(i,j,k)=adv_u(i,j,k)+adv_weno5*dyinv
     
      wcm1=(w(i  ,j  ,k-1)+w(i+1,j  ,k-1))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i+1,j  ,k  ))*0.5d0
      vel00=(wcm1+wcp1)*0.5d0
      dfm5= -u(i,j,k-3)+u(i,j,k-2)
      dfm3= -u(i,j,k-2)+u(i,j,k-1)
      dfm1= -u(i,j,k-1)+u(i,j,k  )
      dfp1= -u(i,j,k  )+u(i,j,k+1)
      dfp3= -u(i,j,k+1)+u(i,j,k+2)
      dfp5= -u(i,j,k+2)+u(i,j,k+3)
      include'weno5.h'
      adv_u(i,j,k)=adv_u(i,j,k)+adv_weno5*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

c
c<v
c
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(dfm5,dfm3,dfm1,dfp1,dfp3,dfp5)
!$OMP$ PRIVATE(vel00)
!$OMP$ PRIVATE(dfm_1,dfm_2,dfm_3)
!$OMP$ PRIVATE(dfp_1,dfp_2,dfp_3)
!$OMP$ PRIVATE(sm_1,sm_2,sm_3)
!$OMP$ PRIVATE(sp_1,sp_2,sp_3)
!$OMP$ PRIVATE(am_1,am_2,am_3)
!$OMP$ PRIVATE(ap_1,ap_2,ap_3)
!$OMP$ PRIVATE(deninv_am,deninv_ap)
!$OMP$ PRIVATE(wgm_1,wgm_3,wgm_2)
!$OMP$ PRIVATE(wgp_1,wgp_2,wgp_3)
!$OMP$ PRIVATE(velab,velm1,velp1,adv_weno5)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
!$OMP$ SHARED(u,v,w,adv_v)
!$OMP$ SHARED(verysmall,d1_12)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      ucm1=(u(i-1,j  ,k  )+u(i-1,j+1,k  ))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i  ,j+1,k  ))*0.5d0
      vel00=(ucm1+ucp1)*0.5d0
      dfm5= -v(i-3,j,k)+v(i-2,j,k)
      dfm3= -v(i-2,j,k)+v(i-1,j,k)
      dfm1= -v(i-1,j,k)+v(i  ,j,k)
      dfp1= -v(i  ,j,k)+v(i+1,j,k)
      dfp3= -v(i+1,j,k)+v(i+2,j,k)
      dfp5= -v(i+2,j,k)+v(i+3,j,k)
      include'weno5.h'
      adv_v(i,j,k)=adv_weno5*dxinv

      vcm1=(v(i  ,j-1,k  )+v(i  ,j  ,k  ))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i  ,j+1,k  ))*0.5d0
      vel00=(vcm1+vcp1)*0.5d0
      dfm5= -v(i,j-3,k)+v(i,j-2,k)
      dfm3= -v(i,j-2,k)+v(i,j-1,k)
      dfm1= -v(i,j-1,k)+v(i,j  ,k)
      dfp1= -v(i,j  ,k)+v(i,j+1,k)
      dfp3= -v(i,j+1,k)+v(i,j+2,k)
      dfp5= -v(i,j+2,k)+v(i,j+3,k)
      include'weno5.h'
      adv_v(i,j,k)=adv_v(i,j,k)+adv_weno5*dyinv

      wcm1=(w(i  ,j  ,k-1)+w(i  ,j+1,k-1))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i  ,j+1,k  ))*0.5d0
      vel00=(wcm1+wcp1)*0.5d0
      dfm5= -v(i,j,k-3)+v(i,j,k-2)
      dfm3= -v(i,j,k-2)+v(i,j,k-1)
      dfm1= -v(i,j,k-1)+v(i,j,k  )
      dfp1= -v(i,j,k  )+v(i,j,k+1)
      dfp3= -v(i,j,k+1)+v(i,j,k+2)
      dfp5= -v(i,j,k+2)+v(i,j,k+3)
      include'weno5.h'
      adv_v(i,j,k)=adv_v(i,j,k)+adv_weno5*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

c
c<w
c
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(dfm5,dfm3,dfm1,dfp1,dfp3,dfp5)
!$OMP$ PRIVATE(vel00)
!$OMP$ PRIVATE(dfm_1,dfm_2,dfm_3)
!$OMP$ PRIVATE(dfp_1,dfp_2,dfp_3)
!$OMP$ PRIVATE(sm_1,sm_2,sm_3)
!$OMP$ PRIVATE(sp_1,sp_2,sp_3)
!$OMP$ PRIVATE(am_1,am_2,am_3)
!$OMP$ PRIVATE(ap_1,ap_2,ap_3)
!$OMP$ PRIVATE(deninv_am,deninv_ap)
!$OMP$ PRIVATE(wgm_1,wgm_3,wgm_2)
!$OMP$ PRIVATE(wgp_1,wgp_2,wgp_3)
!$OMP$ PRIVATE(velab,velm1,velp1,adv_weno5)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
!$OMP$ SHARED(u,v,w,adv_w)
!$OMP$ SHARED(verysmall,d1_12)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      ucm1=(u(i-1,j  ,k  )+u(i-1,j  ,k+1))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i  ,j  ,k+1))*0.5d0
      vel00=(ucm1+ucp1)*0.5d0
      dfm5= -w(i-3,j,k)+w(i-2,j,k)
      dfm3= -w(i-2,j,k)+w(i-1,j,k)
      dfm1= -w(i-1,j,k)+w(i  ,j,k)
      dfp1= -w(i  ,j,k)+w(i+1,j,k)
      dfp3= -w(i+1,j,k)+w(i+2,j,k)
      dfp5= -w(i+2,j,k)+w(i+3,j,k)
      include'weno5.h'
      adv_w(i,j,k)=adv_weno5*dxinv

      vcm1=(v(i  ,j-1,k  )+v(i  ,j-1,k+1))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i  ,j  ,k+1))*0.5d0
      vel00=(vcm1+vcp1)*0.5d0
      dfm5= -w(i,j-3,k)+w(i,j-2,k)
      dfm3= -w(i,j-2,k)+w(i,j-1,k)
      dfm1= -w(i,j-1,k)+w(i,j  ,k)
      dfp1= -w(i,j  ,k)+w(i,j+1,k)
      dfp3= -w(i,j+1,k)+w(i,j+2,k)
      dfp5= -w(i,j+2,k)+w(i,j+3,k)
      include'weno5.h'
      adv_w(i,j,k)=adv_w(i,j,k)+adv_weno5*dyinv

      wcm1=(w(i  ,j  ,k-1)+w(i  ,j  ,k  ))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i  ,j  ,k+1))*0.5d0
      vel00=(wcm1+wcp1)*0.5d0
      dfm5= -w(i,j,k-3)+w(i,j,k-2)
      dfm3= -w(i,j,k-2)+w(i,j,k-1)
      dfm1= -w(i,j,k-1)+w(i,j,k  )
      dfp1= -w(i,j,k  )+w(i,j,k+1)
      dfp3= -w(i,j,k+1)+w(i,j,k+2)
      dfp5= -w(i,j,k+2)+w(i,j,k+3)
      include'weno5.h'
      adv_w(i,j,k)=adv_w(i,j,k)+adv_weno5*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
