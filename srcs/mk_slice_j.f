      subroutine mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk
     & ,j,ip,jp,kp,wrk_slice,qk,qk_slice)
      implicit none
      include 'mpif.h'
      include 'param.h'
      integer ipara,ID,nID(6),ndiv,ni,nj,nk
      integer j,ip,jp,kp
      real*8 wrk_slice(0:ni,0:nk*ndiv)
      real*8        qk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  qk_slice(0:ni,0:nk*ndiv)

      integer nb,ierr
      integer i,k,kst,kk
      
      nb=(ni+1)*(nk*ndiv+1)
      
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk,ndiv)
!$OMP$ SHARED(wrk_slice)
      do k=0,nk*ndiv
      do i=0,ni
      wrk_slice(i,k)=0.0d0
      enddo
      enddo
!$OMP  END PARALLEL DO

      kst=1
      if(nID(Z_MINUS).lt.0)kst=0

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k,kk)
!$OMP$ SHARED(ni,nk,j,kst,ID,ip,jp,kp)
!$OMP$ SHARED(wrk_slice,qk_slice,qk)
      do k=kst,nk
      do i=0,ni
      kk=k+ID*nk
      wrk_slice(i   ,kk)=(
     1       qk(i   ,j   ,k   )
     2 +     qk(i+ip,j   ,k   )
     3 +     qk(i   ,j+jp,k   )
     4 +     qk(i+ip,j+jp,k   )
     5 +     qk(i   ,j   ,k+kp)
     6 +     qk(i+ip,j   ,k+kp)
     7 +     qk(i   ,j+jp,k+kp)
     8 +     qk(i+ip,j+jp,k+kp) )*0.125d0
       qk_slice(i   ,kk)=wrk_slice(i,kk)
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      call mpi_allreduce(
     & wrk_slice,qk_slice,nb,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      return
      end
