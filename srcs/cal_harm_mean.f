ccc
ccc<compute mixture property (harmonic mean) qk from phik using ql and qg
ccc 
      subroutine cal_harm_mean(ni,nj,nk,phik,ql,qg,qk)
      implicit none
      integer ni,nj,nk
      real*8 phik(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 ql,qg
      real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 phil00,phig00
      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,phil00,phig00)
!$OMP$ SHARED(ni,nj,nk,qk,phik,ql,qg)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      phil00=1.0d0-phik(i,j,k)
      phig00=phik(i,j,k)
      qk(i,j,k)=1.0d0/(phil00/ql+phig00/qg)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      return
      end
