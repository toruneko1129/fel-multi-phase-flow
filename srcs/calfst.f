ccc
ccc<compute the surface tension force terms fst_[uvw]k
ccc from the VOF function gradient phix_uk, phiy_vk, phiz_wk
ccc and the curvature rkap_[uvp]k
ccc 

      subroutine calfst(ni,nj,nk
     & ,rhol,rhog,surface_tension
     & ,phix_u,phiy_v,phiz_w,rkap,fst_u,fst_v,fst_w)

      implicit none
      integer ni,nj,nk
      real*8 rhol,rhog,surface_tension
      real*8 phix_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phiy_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phiz_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   rkap(-2:ni+3,-2:nj+3,-2:nk+3,3)
      real*8  fst_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k
      real*8 st_rhominv

      st_rhominv=2.0d0*surface_tension/(rhol+rhog)

ccc
ccc<compute fst_[uvw]
ccc

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(fst_u,fst_v,fst_w,st_rhominv)
!$OMP$ SHARED(rkap,phix_u,phiy_v,phiz_w)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      fst_u(i,j,k)=-st_rhominv*rkap(i,j,k,1)*phix_u(i,j,k)
      fst_v(i,j,k)=-st_rhominv*rkap(i,j,k,2)*phiy_v(i,j,k)
      fst_w(i,j,k)=-st_rhominv*rkap(i,j,k,3)*phiz_w(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

