ccc
ccc<impose boundary conditions on the velocity components uk,vk,wk
ccc
      subroutine bndu(nID,ni,nj,nk,uk,vk,wk)

      implicit none
      include 'param.h'
      integer nID(6)
      integer ni,nj,nk
      real*8    uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    uwall

      integer i,j,k

ccc
ccc<j
ccc
ccc

      uwall = 0.5d0

      if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk)
!$OMP$ SHARED(uk,vk,wk,uwall)
      do k=-2,nk+3
      do i=-2,ni+3
      uk(i,  -2,k)=-uwall*2-uk(i,   3,k)
      uk(i,  -1,k)=-uwall*2-uk(i,   2,k)
      uk(i,   0,k)=-uwall*2-uk(i,   1,k)

      vk(i,  -2,k)=-vk(i,   2,k)
      vk(i,  -1,k)=-vk(i,   1,k)
      vk(i,   0,k)=0.0d0

      wk(i,  -2,k)=-wk(i,   3,k)
      wk(i,  -1,k)=-wk(i,   2,k)
      wk(i,   0,k)=-wk(i,   1,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,vk,wk,uwall)
      do k=-2,nk+3
      do i=-2,ni+3
      uk(i,nj+1,k)=uwall*2-uk(i,nj  ,k)
      uk(i,nj+2,k)=uwall*2-uk(i,nj-1,k)
      uk(i,nj+3,k)=uwall*2-uk(i,nj-2,k)

      vk(i,nj  ,k)=0.0d0
      vk(i,nj+1,k)=-vk(i,nj-1,k)
      vk(i,nj+2,k)=-vk(i,nj-2,k)
      vk(i,nj+3,k)=-vk(i,nj-3,k)

      wk(i,nj+1,k)=-wk(i,nj  ,k)
      wk(i,nj+2,k)=-wk(i,nj-1,k)
      wk(i,nj+3,k)=-wk(i,nj-2,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      return
      end

