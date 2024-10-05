      subroutine mk_slice_k(ipara,ID,nID,ndiv,ni,nj,nk
     & ,k,ip,jp,kp,wrk_slice,qk,qk_slice)
      implicit none
      include 'mpif.h'
      include 'param.h'
      integer ipara,ID,nID(6),ndiv,ni,nj,nk
      integer k,ip,jp,kp
      real*8 wrk_slice(0:ni,0:nj*ndiv)
      real*8        qk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  qk_slice(0:ni,0:nj*ndiv)

      integer nb,ierr
      integer i,j,jst,jj
      
      nb=(ni+1)*(nj*ndiv+1)
      
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j)
!$OMP$ SHARED(ni,nj,ndiv)
!$OMP$ SHARED(wrk_slice)
      do j=0,nj*ndiv
      do i=0,ni
      wrk_slice(i,j)=0.0d0
      enddo
      enddo
!$OMP  END PARALLEL DO

      jst=1
      if(nID(Y_MINUS).lt.0)jst=0

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,jj)
!$OMP$ SHARED(ni,nj,k,jst,ID,ip,jp,kp)
!$OMP$ SHARED(wrk_slice,qk_slice,qk)
      do j=jst,nj
      do i=0,ni
      jj=j+ID*nj
      wrk_slice(i   ,jj)=(
     1       qk(i   ,j   ,k   )
     2 +     qk(i+ip,j   ,k   )
     3 +     qk(i   ,j+jp,k   )
     4 +     qk(i+ip,j+jp,k   )
     5 +     qk(i   ,j   ,k+kp)
     6 +     qk(i+ip,j   ,k+kp)
     7 +     qk(i   ,j+jp,k+kp)
     8 +     qk(i+ip,j+jp,k+kp) )*0.125d0
       qk_slice(i   ,jj)=wrk_slice(i,jj)
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      call mpi_allreduce(
     & wrk_slice,qk_slice,nb,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      return
      end
