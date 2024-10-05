ccc
ccc<impose the Neumann boundary condition on qk
ccc
      subroutine bnd_neumann(nID,ni,nj,nk,qk)

      implicit none
      include 'param.h'
      integer nID(6),ni,nj,nk
      real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

      if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i,  -2,k)=qk(i,   3,k)
      qk(i,  -1,k)=qk(i,   2,k)
      qk(i,   0,k)=qk(i,   1,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nj,nk,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i,nj+1,k)=qk(i,nj  ,k)
      qk(i,nj+2,k)=qk(i,nj-1,k)
      qk(i,nj+3,k)=qk(i,nj-2,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      return
      end
