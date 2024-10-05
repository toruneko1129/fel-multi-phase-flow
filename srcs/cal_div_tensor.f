      subroutine cal_div_tensor(ni,nj,nk,dxinv,dyinv,dzinv
     & ,tauk,uk,vk,wk)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8 tauk(-2:ni+3,-2:nj+3,-2:nk+3,6)
      real*8   uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   wk(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,vk,wk)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      uk(i,j,k)=0.0d0
      vk(i,j,k)=0.0d0
      wk(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,tauk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+3
      do j=-1,nj+3
      do i=-2,ni+2
      uk(i,j,k)=
     1 +(-tauk(i,j  ,k  ,1)+tauk(i+1,j,k,1))*dxinv
     2 +(-tauk(i,j-1,k  ,4)+tauk(i  ,j,k,4))*dyinv
     3 +(-tauk(i,j  ,k-1,5)+tauk(i  ,j,k,5))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(vk,tauk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-1,nk+3
      do j=-2,nj+2
      do i=-1,ni+3
      vk(i,j,k)=
     1 +(-tauk(i-1,j,k  ,4)+tauk(i,j  ,k,4))*dxinv
     2 +(-tauk(i  ,j,k  ,2)+tauk(i,j+1,k,2))*dyinv
     3 +(-tauk(i  ,j,k-1,6)+tauk(i,j  ,k,6))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(wk,tauk)
!$OMP$ SHARED(dxinv,dyinv,dzinv)
      do k=-2,nk+2
      do j=-1,nj+3
      do i=-1,ni+3
      wk(i,j,k)=
     1 +(-tauk(i-1,j  ,k,5)+tauk(i,j,k  ,5))*dxinv
     2 +(-tauk(i  ,j-1,k,6)+tauk(i,j,k  ,6))*dyinv
     3 +(-tauk(i  ,j  ,k,3)+tauk(i,j,k+1,3))*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

