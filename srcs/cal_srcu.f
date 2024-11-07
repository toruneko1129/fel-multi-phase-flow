      subroutine cal_srcu(nID,ni,nj,nk
     & ,rho_av,dt,grvb,grvp
     & ,u,v,w,rho,rhon
     & ,adv_u ,adv_v ,adv_w
     & ,adv_uo,adv_vo,adv_wo
     & ,prs_u ,prs_v ,prs_w 
     & ,vis_u ,vis_v ,vis_w 
     & ,fst_u ,fst_v ,fst_w
     & ,fst_un,fst_vn,fst_wn
     & ,src_u ,src_v ,src_w )

      implicit none
      include 'param.h'
      integer nID(6)
      integer ni,nj,nk
      real*8 rho_av,dt,grvb,grvp
      real*8       u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8       v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8       w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     rho(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    rhon(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   adv_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   adv_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   adv_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_uo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_vo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  adv_wo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   prs_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   prs_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   prs_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vis_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vis_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vis_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   fst_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   fst_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   fst_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_un(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_vn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_wn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k,jen
      real*8 dtinv,rho00

      jen=nj
      if(nID(Y_PLUS).lt.0)jen=nj-1

      dtinv=1.0d0/dt


!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rho00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rho,rhon,src_u,u,dtinv,adv_u,adv_uo)
!$OMP$ SHARED(prs_u,vis_u,fst_u,fst_un,grvb,rho_av)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      rho00=(
     1   rho(i,j,k)+ rho(i+1,j,k)
     2 +rhon(i,j,k)+rhon(i+1,j,k))*0.25d0
         src_u(i,j,k)=
     1 +     u(i,j,k)*dtinv
     2 - adv_u(i,j,k)*1.5d0
     3 +adv_uo(i,j,k)*0.5d0
     4 - prs_u(i,j,k)/rho00
     5 + vis_u(i,j,k)/rho00*0.5d0
     6 + fst_u(i,j,k)*0.5d0
     7 +fst_un(i,j,k)*0.5d0
     8 + grvb
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!+ grvb*(rho_av/rho00-1.0d0)

!6 + fst_u(i,j,k)*0.5d0
!7 +fst_un(i,j,k)*0.5d0

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rho00)
!$OMP$ SHARED(ni,jen,nk)
!$OMP$ SHARED(rho,rhon,src_v,v,dtinv,adv_v,adv_vo)
!$OMP$ SHARED(prs_v,vis_v,fst_v,fst_vn,grvp,rho_av)
      do k=1,nk
      do j=1,jen
      do i=1,ni
      rho00=(
     1   rho(i,j,k)+ rho(i,j+1,k)
     2 +rhon(i,j,k)+rhon(i,j+1,k))*0.25d0
         src_v(i,j,k)=
     1 +     v(i,j,k)*dtinv
     2 - adv_v(i,j,k)*1.5d0
     3 +adv_vo(i,j,k)*0.5d0
     4 - prs_v(i,j,k)/rho00
     5 + vis_v(i,j,k)/rho00*0.5d0
     6 + fst_v(i,j,k)*0.5d0
     7 +fst_vn(i,j,k)*0.5d0
     8 + grvp
       enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!+ grvp*(rho_av/rho00-1.0d0)
!6 + fst_v(i,j,k)*0.5d0
!7 +fst_vn(i,j,k)*0.5d0

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rho00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rho,rhon,src_w,w,dtinv,adv_w,adv_wo)
!$OMP$ SHARED(prs_w,vis_w,fst_w,fst_wn)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      rho00=(
     1   rho(i,j,k)+ rho(i,j,k+1)
     2 +rhon(i,j,k)+rhon(i,j,k+1))*0.25d0
         src_w(i,j,k)=
     1 +     w(i,j,k)*dtinv
     2 - adv_w(i,j,k)*1.5d0
     3 +adv_wo(i,j,k)*0.5d0
     4 - prs_w(i,j,k)/rho00
     5 + vis_w(i,j,k)/rho00*0.5d0
     6 + fst_w(i,j,k)*0.5d0
     7 +fst_wn(i,j,k)*0.5d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
!6 + fst_w(i,j,k)*0.5d0
!7 +fst_wn(i,j,k)*0.5d0
      return
      end
