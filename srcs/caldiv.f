      subroutine caldiv(ipara,ndiv,ni,nj,nk,dxinv,dyinv,dzinv,dt
     & ,uk,vk,wk,div,err_div)

      implicit none
      include 'mpif.h'
      integer ipara,ndiv
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv,dt
      real*8  uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 div(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 err_div

      integer ierr
      real*8 sendbuf(1),recvbuf(1)
      integer i,j,k
      real*8 dtinv

cc
cc<compute divergence divided by dt
cc

      dtinv=1.0d0/dt
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(div,uk,vk,wk,dxinv,dyinv,dzinv,dtinv)
      do k=-1,nk+3
      do j=-1,nj+3
      do i=-1,ni+3
      div(i,j,k)=(
     1  (-uk(i-1,j  ,k  )+uk(i,j,k))*dxinv
     2 +(-vk(i  ,j-1,k  )+vk(i,j,k))*dyinv
     2 +(-wk(i  ,j  ,k-1)+wk(i,j,k))*dzinv
     & )*dtinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      err_div=0.0d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(div)
!$OMP$ REDUCTION(+:err_div)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      err_div=err_div+div(i,j,k)**2
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      sendbuf(1) = err_div
      recvbuf(1) = err_div
      call mpi_allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      err_div = recvbuf(1)
      endif

      err_div=sqrt(err_div/dble(ni*nj*nk*ndiv))+1.0d-99

      return
      end
