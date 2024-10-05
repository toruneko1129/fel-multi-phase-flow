        subroutine sum_fst(ni,nj,nk
     & ,fst_u,fst_v,fst_w
     & ,sum_fst_u,sum_fst_v,sum_fst_w)

      implicit none
      integer ni,nj,nk
      real*8  fst_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  fst_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  sum_fst_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  sum_fst_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  sum_fst_w(-2:ni+3,-2:nj+3,-2:nk+3)
      integer i,j,k

ccc
ccc<summation fst_[uvw] fst_[uvw]n
ccc

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(fst_u,fst_v,fst_w)
!$OMP$ SHARED(sum_fst_u,sum_fst_v,sum_fst_w)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      sum_fst_u(i,j,k)=sum_fst_u(i,j,k)+fst_u(i,j,k)
      sum_fst_v(i,j,k)=sum_fst_v(i,j,k)+fst_v(i,j,k)
      sum_fst_w(i,j,k)=sum_fst_w(i,j,k)+fst_w(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
