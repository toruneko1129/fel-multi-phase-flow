ccc
ccc<compute expansion coefficients for pressure equation using fft
ccc 

      subroutine cal_coef_solp_fft(nj,as_p,an_p
     & ,as_p_fft,an_p_fft,ap_p_fft)

      implicit none
      integer nj
      real*8     as_p(-2:nj+3)
      real*8     an_p(-2:nj+3)
      real*8 as_p_fft(-2:nj+3)
      real*8 an_p_fft(-2:nj+3)
      real*8 ap_p_fft(-2:nj+3)

      integer j

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(j)
!$OMP$ SHARED(nj)
!$OMP$ SHARED(as_p,an_p)
!$OMP$ SHARED(as_p_fft,an_p_fft,ap_p_fft)
      do j=1,nj
      as_p_fft(j)=as_p(j)
      an_p_fft(j)=an_p(j)
      ap_p_fft(j)=as_p_fft(j)+an_p_fft(j)
      enddo
!$OMP  END PARALLEL DO


c< for applying the Neumann condition on the wall

      ap_p_fft(1)=ap_p_fft(1)-as_p_fft(1)
      as_p_fft(1)=0.0d0

      ap_p_fft(nj)=ap_p_fft(nj)-an_p_fft(nj)
      an_p_fft(nj)=0.0d0

      return
      end
