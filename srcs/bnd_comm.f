      subroutine bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,qk)

      implicit none
      include 'mpif.h'
      include 'param.h'
      integer ipara,nID(6)
      integer ni,nj,nk
      integer  key(2,2)
      real*8 recvjb(-2:ni+3,-2:nk+3,3,2)
      real*8 sendjb(-2:ni+3,-2:nk+3,3,2)
      real*8     qk(-2:ni+3,-2:nj+3,-2:nk+3)

      integer status(MPI_STATUS_SIZE)
      integer ierr
      integer i,k
      integer nb

      if(ipara.eq.1)then
      nb=(ni+6)*(nk+6)*3
      key(1,1)=MPI_REQUEST_NULL
      key(2,1)=MPI_REQUEST_NULL
      key(1,2)=MPI_REQUEST_NULL
      key(2,2)=MPI_REQUEST_NULL

c<j
      if(nID(Y_MINUS).ge.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk)
!$OMP$ SHARED(sendjb,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      sendjb(i,k,1,1)=qk(i,   1,k)
      sendjb(i,k,2,1)=qk(i,   2,k)
      sendjb(i,k,3,1)=qk(i,   3,k)
      enddo
      enddo
!$OMP END PARALLEL DO
      call mpi_isend(sendjb(-2,-2,1,1),nb,MPI_REAL8,
     &  nID(Y_MINUS),0,MPI_COMM_WORLD,key(1,1),ierr)
      endif

      if(nID(Y_PLUS ).ge.0)then
      call mpi_irecv(recvjb(-2,-2,1,1),nb,MPI_REAL8,
     &  nID(Y_PLUS ),0,MPI_COMM_WORLD,key(1,2),ierr)

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(sendjb,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      sendjb(i,k,1,2)=qk(i,nj-2,k)
      sendjb(i,k,2,2)=qk(i,nj-1,k)
      sendjb(i,k,3,2)=qk(i,nj  ,k)
      enddo
      enddo
!$OMP END PARALLEL DO
      call mpi_isend(sendjb(-2,-2,1,2),nb,MPI_REAL8,
     &  nID(Y_PLUS ),0,MPI_COMM_WORLD,key(2,1),ierr)
      endif

      if(nID(Y_MINUS).ge.0)then
      call mpi_irecv(recvjb(-2,-2,1,2),nb,MPI_REAL8,
     &  nID(Y_MINUS),0,MPI_COMM_WORLD,key(2,2),ierr)
      endif

      if(           key(1,1).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(1,1),status,ierr)
      endif

      if(           key(1,2).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(1,2),status,ierr)
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(recvjb,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i,nj+1,k)=recvjb(i,k,1,1)
      qk(i,nj+2,k)=recvjb(i,k,2,1)
      qk(i,nj+3,k)=recvjb(i,k,3,1)
      enddo
      enddo
!$OMP END PARALLEL DO
      endif

      if(           key(2,1).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(2,1),status,ierr)
      endif

      if(           key(2,2).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(2,2),status,ierr)
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk)
!$OMP$ SHARED(recvjb,qk)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i,  -2,k)=recvjb(i,k,1,2)
      qk(i,  -1,k)=recvjb(i,k,2,2)
      qk(i,   0,k)=recvjb(i,k,3,2)
      enddo
      enddo
!$OMP END PARALLEL DO
      endif

      endif


      return
      end

