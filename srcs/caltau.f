ccc
ccc<compute the viscous stress components: tau[ij]k
ccc from the strain rate components s[ij]k and the viscosity rmu
ccc

      subroutine caltau(ni,nj,nk,rmu,sk,tauk)
      implicit none
      integer ni,nj,nk
      real*8 rmu
      real*8   sk(-2:ni+3,-2:nj+3,-2:nk+3,6)
      real*8 tauk(-2:ni+3,-2:nj+3,-2:nk+3,6)


      integer i,j,k,l

      do l=1,6
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(l,ni,nj,nk)
!$OMP$ SHARED(tauk,sk)
!$OMP$ SHARED(rmu)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      tauk(i,j,k,l)=sk(i,j,k,l)*rmu*2.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

      return
      end
