      subroutine mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,qk,qk_all)
      implicit none
      include 'mpif.h'
      include 'param.h'
      integer ipara,ID,nID(6),ndiv,ni,nj,nk
      real*8 wrk_all(-2:ni+3,-2:nj*ndiv+3,-2:nk+3)
      real*8      qk(-2:ni+3,-2:nj     +3,-2:nk+3)
      real*8  qk_all(-2:ni+3,-2:nj*ndiv+3,-2:nk+3)

      integer nb,ierr
      integer i,j,k,jst,jen,jj

      nb=(ni+6)*(nj*ndiv+6)*(nk+6)

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,ndiv)
!$OMP$ SHARED(wrk_all)
      do k=-2,nk     +3
      do j=-2,nj*ndiv+3
      do i=-2,ni     +3
      wrk_all(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      jst=1
      if(nID(Y_MINUS).lt.0)jst=-2
      jen=nj
      if(nID(Y_PLUS ).lt.0)jen=nj+3

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,jj)
!$OMP$ SHARED(ni,nj,nk,jst,jen,ID)
!$OMP$ SHARED(wrk_all,qk_all,qk)
      do k=-2,nk+3
      do j=jst,jen
      do i=-2,ni+3
      jj=j+ID*nj
      wrk_all(i,jj,k)=qk(i,j,k)
       qk_all(i,jj,k)=qk(i,j,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      call mpi_allreduce(
     & wrk_all,qk_all,nb,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      return
      end
