      subroutine cal_av(ipara,ndiv,ni,nj,nk,q,av)
      implicit none
      include 'mpif.h'
      integer ipara,ndiv,ni,nj,nk
      real*8 q(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av

      integer ierr
      real*8 sendbuf(1),recvbuf(1)
      integer i,j,k
      real*8 sum

      sum=0.0d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(q)
!$OMP$ REDUCTION(+:sum)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      sum=sum+q(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      sendbuf(1) = sum
      recvbuf(1) = sum
      call mpi_allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      sum = recvbuf(1)
      endif

      av=sum/dble(ni*nj*nk*ndiv)

      return
      end
