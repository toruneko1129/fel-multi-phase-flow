      subroutine calrkap(ni,nj,nk,dxinv,dyinv,dzinv
     & ,rn_pp,rn_12,rn_13,rn_23,rkap)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8 rn_pp(-2:ni+3,-2:nj+3,-2:nk+3,3)
      real*8 rn_12(-2:ni+3,-2:nj+3,-2:nk+3,3)
      real*8 rn_13(-2:ni+3,-2:nj+3,-2:nk+3,3)
      real*8 rn_23(-2:ni+3,-2:nj+3,-2:nk+3,3)
      real*8  rkap(-2:ni+3,-2:nj+3,-2:nk+3,3)

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,rkap)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      rkap(i,j,k,1)=0.0d0
      rkap(i,j,k,2)=0.0d0
      rkap(i,j,k,3)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<compute curvature rkap(i,j,k,1) at the definition point of u
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rkap,rn_pp,rn_12,rn_13)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+3
      do j=-1,nj+3
      do i=-2,ni+2
      rkap(i,j,k,1)=
     1 +(-rn_pp(i,j  ,k  ,1)+rn_pp(i+1,j,k,1))*dxinv
     2 +(-rn_12(i,j-1,k  ,2)+rn_12(i  ,j,k,2))*dyinv
     3 +(-rn_13(i,j  ,k-1,3)+rn_13(i  ,j,k,3))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<compute curvature rkap(i,j,k,2) at the definition point of v
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rkap,rn_12,rn_pp,rn_23)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+3
      do j=-2,nj+2
      do i=-1,ni+3
      rkap(i,j,k,2)=
     1 +(-rn_12(i-1,j,k  ,1)+rn_12(i,j  ,k,1))*dxinv
     2 +(-rn_pp(i  ,j,k  ,2)+rn_pp(i,j+1,k,2))*dyinv
     3 +(-rn_23(i  ,j,k-1,3)+rn_23(i,j  ,k,3))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<compute curvature rkap(i,j,k,3) at the definition point of w
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rkap,rn_13,rn_23,rn_pp)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-2,nk+2
      do j=-1,nj+3
      do i=-1,ni+3
      rkap(i,j,k,3)=
     1 +(-rn_13(i-1,j  ,k,1)+rn_13(i,j,k  ,1))*dxinv
     2 +(-rn_23(i  ,j-1,k,2)+rn_23(i,j,k  ,2))*dyinv
     3 +(-rn_pp(i  ,j  ,k,3)+rn_pp(i,j,k+1,3))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

