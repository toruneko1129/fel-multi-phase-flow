ccc
ccc<compute the advection terms adv_uk, adv_vk, adv_wk
ccc from the velocity vector components uk, vk, wk
ccc 

      subroutine cal_advu(ni,nj,nk,dxinv,dyinv,dzinv
     & ,u,v,w,adv_u,adv_v,adv_w)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8     u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 adv_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 adv_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 adv_w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k
      real*8 ucm1,ucp1,vcm1,vcp1,wcm1,wcp1
      real*8 gxm1,gxp1,gym1,gyp1,gzm1,gzp1

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

ccc
ccc<compute adv_u
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(gxm1,gxp1,gym1,gyp1,gzm1,gzp1)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,v,w,adv_u)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      ucm1=(u(i-1,j  ,k  )+u(i  ,j  ,k  ))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i+1,j  ,k  ))*0.5d0
      vcm1=(v(i  ,j-1,k  )+v(i+1,j-1,k  ))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i+1,j  ,k  ))*0.5d0
      wcm1=(w(i  ,j  ,k-1)+w(i+1,j  ,k-1))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i+1,j  ,k  ))*0.5d0

      gxm1=-u(i-1,j  ,k  )+u(i  ,j  ,k  )
      gxp1=-u(i  ,j  ,k  )+u(i+1,j  ,k  )
      gym1=-u(i  ,j-1,k  )+u(i  ,j  ,k  )
      gyp1=-u(i  ,j  ,k  )+u(i  ,j+1,k  )
      gzm1=-u(i  ,j  ,k-1)+u(i  ,j  ,k  )
      gzp1=-u(i  ,j  ,k  )+u(i  ,j  ,k+1)

      adv_u(i,j,k)=
     1 +(ucm1*gxm1+ucp1*gxp1)*0.5d0*dxinv
     2 +(vcm1*gym1+vcp1*gyp1)*0.5d0*dyinv
     3 +(wcm1*gzm1+wcp1*gzp1)*0.5d0*dzinv

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc   
ccc<compute adv_v
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(gxm1,gxp1,gym1,gyp1,gzm1,gzp1)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,v,w,adv_v)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      ucm1=(u(i-1,j  ,k  )+u(i-1,j+1,k  ))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i  ,j+1,k  ))*0.5d0
      vcm1=(v(i  ,j-1,k  )+v(i  ,j  ,k  ))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i  ,j+1,k  ))*0.5d0
      wcm1=(w(i  ,j  ,k-1)+w(i  ,j+1,k-1))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i  ,j+1,k  ))*0.5d0

      gxm1=-v(i-1,j  ,k  )+v(i  ,j  ,k  )
      gxp1=-v(i  ,j  ,k  )+v(i+1,j  ,k  )
      gym1=-v(i  ,j-1,k  )+v(i  ,j  ,k  )
      gyp1=-v(i  ,j  ,k  )+v(i  ,j+1,k  )
      gzm1=-v(i  ,j  ,k-1)+v(i  ,j  ,k  )
      gzp1=-v(i  ,j  ,k  )+v(i  ,j  ,k+1)

      adv_v(i,j,k)=
     1 +(ucm1*gxm1+ucp1*gxp1)*0.5d0*dxinv
     2 +(vcm1*gym1+vcp1*gyp1)*0.5d0*dyinv
     3 +(wcm1*gzm1+wcp1*gzp1)*0.5d0*dzinv

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<compute adv_w
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(ucm1,ucp1,vcm1,vcp1,wcm1,wcp1)
!$OMP$ PRIVATE(gxm1,gxp1,gym1,gyp1,gzm1,gzp1)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,v,w,adv_w)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      ucm1=(u(i-1,j  ,k  )+u(i-1,j  ,k+1))*0.5d0
      ucp1=(u(i  ,j  ,k  )+u(i  ,j  ,k+1))*0.5d0
      vcm1=(v(i  ,j-1,k  )+v(i  ,j-1,k+1))*0.5d0
      vcp1=(v(i  ,j  ,k  )+v(i  ,j  ,k+1))*0.5d0
      wcm1=(w(i  ,j  ,k-1)+w(i  ,j  ,k  ))*0.5d0
      wcp1=(w(i  ,j  ,k  )+w(i  ,j  ,k+1))*0.5d0

      gxm1=-w(i-1,j  ,k  )+w(i  ,j  ,k  )
      gxp1=-w(i  ,j  ,k  )+w(i+1,j  ,k  )
      gym1=-w(i  ,j-1,k  )+w(i  ,j  ,k  )
      gyp1=-w(i  ,j  ,k  )+w(i  ,j+1,k  )
      gzm1=-w(i  ,j  ,k-1)+w(i  ,j  ,k  )
      gzp1=-w(i  ,j  ,k  )+w(i  ,j  ,k+1)

      adv_w(i,j,k)=
     1 +(ucm1*gxm1+ucp1*gxp1)*0.5d0*dxinv
     2 +(vcm1*gym1+vcp1*gyp1)*0.5d0*dyinv
     3 +(wcm1*gzm1+wcp1*gzp1)*0.5d0*dzinv

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

