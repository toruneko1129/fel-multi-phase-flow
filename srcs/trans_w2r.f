      subroutine trans_w2r(ipara,ID,ndiv,svall,key,sendbuf,recvbuf,w,r)
      implicit none
      include 'mpif.h'

      integer ipara,ID,ndiv,svall(3)
      integer  key(0:ndiv-1,0:ndiv-1,2)
      real*8 sendbuf(-svall(1)/2:svall(1)/2-1,svall(2)/ndiv
     & ,svall(3)/ndiv,0:ndiv-1)
      real*8 recvbuf(-svall(1)/2:svall(1)/2-1,svall(2)/ndiv
     & ,svall(3)/ndiv,0:ndiv-1)
      real*8       w(-svall(1)/2:svall(1)/2-1,svall(2)     
     & , svall(3)/ndiv    )
      real*8       r(-svall(1)/2:svall(1)/2-1,svall(2)/ndiv
     & ,-svall(3)/2:svall(3)/2-1)

      integer status(MPI_STATUS_SIZE)
      integer nb,ierr
      integer i,j,k,l,jj,kk

      nb=svall(1)*svall(2)*svall(3)/ndiv/ndiv

c<for serial computing
      if(ipara.eq.0)then

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,kk)
!$OMP$ SHARED(svall,r,w)
      do k=          1,svall(3)
      do j=          1,svall(2)
      do i=-svall(1)/2,svall(1)/2-1
      kk=k-svall(3)/2-1
      r(i,j,kk)=w(i,j,k)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

c<for MPI parallel computing
      else

c<send and receive data
      do l=0,ndiv-1
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,jj)
!$OMP$ SHARED(l,ndiv,svall)
!$OMP$ SHARED(sendbuf,recvbuf,w)
      do k=          1,svall(3)/ndiv
      do j=          1,svall(2)/ndiv
      do i=-svall(1)/2,svall(1)/2-1
      jj=j+l*svall(2)/ndiv
      sendbuf(i,j,k,l)=w(i,jj,k)
      recvbuf(i,j,k,l)=w(i,jj,k)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo

      do l=0,ndiv-1
      key(ID,l,1)=MPI_REQUEST_NULL
      key(ID,l,2)=MPI_REQUEST_NULL
      enddo

      do l=0,ndiv-1
      if(l.ne.ID)then
      call mpi_isend(
     &  sendbuf(-svall(1)/2,1,1,l),nb,MPI_REAL8
     & ,l,0,MPI_COMM_WORLD,key(ID,l,1),ierr)
      call mpi_irecv(
     &  recvbuf(-svall(1)/2,1,1,l),nb,MPI_REAL8
     & ,l,0,MPI_COMM_WORLD,key(ID,l,2),ierr)
      if(           key(ID,l,1).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(ID,l,1),status,ierr)
      endif
      if(           key(ID,l,2).ne.MPI_REQUEST_NULL)then
      call mpi_wait(key(ID,l,2),status,ierr)
      endif
      endif
      enddo

      do l=0,ndiv-1
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,kk)
!$OMP$ SHARED(l,ndiv,svall)
!$OMP$ SHARED(r,recvbuf)
      do k=          1,svall(3)/ndiv
      do j=          1,svall(2)/ndiv
      do i=-svall(1)/2,svall(1)/2-1
      kk=k+l*svall(3)/ndiv-svall(3)/2-1
      r(i,j,kk)=recvbuf(i,j,k,l)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo

      endif

      return
      end

      


     
