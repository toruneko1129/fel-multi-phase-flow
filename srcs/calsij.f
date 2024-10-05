cccc
ccc<compute strain rate tensor sk
ccc

      subroutine calsij(ni,nj,nk,dxinv,dyinv,dzinv
     & ,uk,vk,wk,sk)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8   uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   sk(-2:ni+3,-2:nj+3,-2:nk+3,6)

      integer i,j,k

ccc
ccc<s11k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,uk,dxinv)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-1,ni+3
      sk(i,j,k,1)=(-uk(i-1,j,k)+uk(i,j,k))*dxinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<s22k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,vk,dyinv)
      do k=-2,nk+3
      do j=-1,nj+3
      do i=-2,ni+3
      sk(i,j,k,2)=(-vk(i,j-1,k)+vk(i,j,k))*dyinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<s33k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,wk,dzinv)
      do k=-1,nk+3
      do j=-1,nj+3
      do i=-1,ni+3
      sk(i,j,k,3)=(-wk(i,j,k-1)+wk(i,j,k))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<s12k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,uk,vk,dxinv,dyinv)
      do k=-2,nk+3
      do j=-2,nj+2
      do i=-2,ni+2
      sk(i,j,k,4)=(
     1 +(-uk(i,j,k)+uk(i  ,j+1,k))*dyinv
     2 +(-vk(i,j,k)+vk(i+1,j  ,k))*dxinv
     & )*0.5d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<s13k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,uk,wk,dxinv,dzinv)
      do k=-2,nk+2
      do j=-2,nj+3
      do i=-2,ni+2
      sk(i,j,k,5)=(
     1 +(-uk(i,j,k)+uk(i  ,j,k+1))*dzinv
     2 +(-wk(i,j,k)+wk(i+1,j,k  ))*dxinv
     & )*0.5d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<s23k
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sk,vk,wk,dyinv,dzinv)
      do k=-2,nk+2
      do j=-2,nj+2
      do i=-2,ni+3
      sk(i,j,k,6)=(
     1  (-vk(i,j,k)+vk(i,j  ,k+1))*dzinv
     2 +(-wk(i,j,k)+wk(i,j+1,k  ))*dyinv
     & )*0.5d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO


      return
      end
