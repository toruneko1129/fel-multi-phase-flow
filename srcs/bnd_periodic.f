ccc
ccc<impose periodic boundary condition on qk
ccc
      subroutine bnd_periodic(ni,nj,nk,qk)
      implicit none
      integer ni,nj,nk
      real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(j,k)
!$OMP$ SHARED(ni,nj,nk,qk)
      do k=-2,nk+3
      do j=-2,nj+3
      qk(  -2,j,k)=qk(ni-2,j,k)
      qk(  -1,j,k)=qk(ni-1,j,k)
      qk(   0,j,k)=qk(ni  ,j,k)
      qk(ni+1,j,k)=qk(   1,j,k)
      qk(ni+2,j,k)=qk(   2,j,k)
      qk(ni+3,j,k)=qk(   3,j,k)
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j)
!$OMP$ SHARED(ni,nj,nk,qk)
      do j=-2,nj+3
      do i=-2,ni+3
      qk(i,j,  -2)=qk(i,j,nk-2)
      qk(i,j,  -1)=qk(i,j,nk-1)
      qk(i,j,   0)=qk(i,j,nk  )
      qk(i,j,nk+1)=qk(i,j,   1)
      qk(i,j,nk+2)=qk(i,j,   2)
      qk(i,j,nk+3)=qk(i,j,   3)
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

