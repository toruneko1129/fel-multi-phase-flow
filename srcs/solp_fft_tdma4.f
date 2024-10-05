ccc
ccc solve pressure equation using FFT and TDMA
ccc
      subroutine solp_fft_tdma4(ipara,ID,ndiv,ni,nj,nk,nstep,imon_t
     & ,rho,dxinv,dyinv,dzinv,div,dp)

      implicit none
      include 'mpif.h'
      integer ipara,ID,ndiv,ni,nj,nk,nstep,imon_t
      real*8  rho,dxinv,dyinv,dzinv
      real*8      div(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8       dp(-2:ni+3,-2:nj+3,-2:nk+3)
      
      integer ierr
      real*8 sendbuf(2),recvbuf(2)
      integer i,j,k
      real*8 aw,ae,as,an,ab,at,ap
      real*8 err0,err,ddp

      aw=dxinv**2/rho
      ae=dxinv**2/rho
      as=dyinv**2/rho
      an=dyinv**2/rho
      ab=dzinv**2/rho
      at=dzinv**2/rho
      ap=aw+ae+as+an+ab+at

      if(mod(nstep,imon_t).eq.0)then

      err0=0.0d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(div)
!$OMP$ REDUCTION(+:err0)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      err0=err0+div(i,j,k)**2
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      err=0.0d0
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,ddp)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(dp,div)
!$OMP$ SHARED(aw,ae,as,an,ab,at,ap)
!$OMP$ REDUCTION(+:err)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      ddp=      
     1 + dp(i  ,j  ,k-1)*ab
     1 + dp(i  ,j-1,k  )*as
     2 + dp(i-1,j  ,k  )*aw
     3 - dp(i  ,j  ,k  )*ap
     4 + dp(i+1,j  ,k  )*ae
     5 + dp(i  ,j+1,k  )*an
     1 + dp(i  ,j  ,k+1)*at
     6 -div(i  ,j  ,k  )
      err=err+ddp**2
c      write(100,'(20e20.10)')dble(i),dble(j),dble(k),ddp
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      sendbuf(1) = err0
      recvbuf(1) = err0
      sendbuf(2) = err
      recvbuf(2) = err
      call mpi_allreduce(sendbuf,recvbuf,2,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      err0 = recvbuf(1)
      err  = recvbuf(2)
      endif
      err0=sqrt(err0/dble(ni*nj*nk*ndiv))+1.0d-99
      err =sqrt(err /dble(ni*nj*nk*ndiv))+1.0d-99

      if(ID.eq.0)then
      write(*,'("Err_solp ",1i10,20e20.10)')
     2  nstep
     3 ,err/err0
     4 ,err
     5 ,err0
      endif

      endif

      return
      end

c -------------------------------------------------------------

