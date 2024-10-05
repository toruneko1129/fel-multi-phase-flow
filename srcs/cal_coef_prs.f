ccc
ccc<compute expansion coefficients for pressure equation 
ccc 

      subroutine cal_coef_prs(nj,rho,dxinv,dyinv,dzinv
     & ,aw_pk,ae_pk,as_pk,an_pk,ab_pk,at_pk,ap_pk)

      implicit none
      integer nj
      real*8 rho,dxinv,dyinv,dzinv
      real*8 aw_pk(-2:nj+3)
      real*8 ae_pk(-2:nj+3)
      real*8 as_pk(-2:nj+3)
      real*8 an_pk(-2:nj+3)
      real*8 ab_pk(-2:nj+3)
      real*8 at_pk(-2:nj+3)
      real*8 ap_pk(-2:nj+3)

      integer j

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(j)
!$OMP$ SHARED(nj)
!$OMP$ SHARED(rho,dxinv,dyinv,dzinv)
!$OMP$ SHARED(aw_pk,ae_pk,as_pk,an_pk,ab_pk,at_pk,ap_pk)
      do j=1,nj

      aw_pk(j)=dxinv**2/rho
      ae_pk(j)=dxinv**2/rho
      as_pk(j)=dyinv**2/rho
      an_pk(j)=dyinv**2/rho
      ab_pk(j)=dzinv**2/rho
      at_pk(j)=dzinv**2/rho

        ap_pk(j)=
     1 +aw_pk(j)
     2 +ae_pk(j)
     3 +as_pk(j)
     4 +an_pk(j)
     5 +ab_pk(j)
     6 +at_pk(j)
      enddo
!$OMP  END PARALLEL DO


c< for applying the Neumann condition on the walls

      ap_pk(1)=ap_pk(1)-as_pk(1)
      as_pk(1)=0.0d0

      ap_pk(nj)=ap_pk(nj)-an_pk(nj)
      an_pk(nj)=0.0d0

      return
      end
