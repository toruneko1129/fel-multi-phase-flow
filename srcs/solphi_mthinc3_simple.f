      subroutine solphi_mthinc3(ipara,ni,nj,nk
     & ,dxinv,dyinv,dzinv
     & ,bet_mthinc
     & ,phix,phiy,phiz
     & ,phi,phin)

      implicit none
      include 'mpif.h'
      integer ipara
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8 bet_mthinc
      real*8   phix(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phiy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phiz(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phin(-2:ni+3,-2:nj+3,-2:nk+3)

      integer ierr
      real*8 sendbuf(4),recvbuf(4)
      integer i,j,k
      real*8 dx,dy,dz
      real*8 bet,dcellinv,eps_mthinc,eps_trunc
      real*8 rnx00,rny00,rnz00,q00
      real*8 sum1,sum2,den1,den2,offset

      dx=1.0d0/dxinv
      dy=1.0d0/dyinv
      dz=1.0d0/dzinv
      bet=bet_mthinc
      dcellinv=1.0d0/sqrt((dx**2+dy**2+dz**2)/3.0d0)
      eps_mthinc=bet/(2.0d0*cosh(-bet*1.5d0)**2)
      eps_trunc=0.5d0*(tanh(-bet*2.5d0)+1.0d0)

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(rnx00,rny00,rnz00,q00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(phix,phiy,phiz)
!$OMP$ SHARED(eps_mthinc,dcellinv,phin,eps_trunc)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1

         rnx00=(
     1 +phix(i-1,j-1,k-1)
     2 +phix(i  ,j-1,k-1)
     3 +phix(i-1,j  ,k-1)
     4 +phix(i  ,j  ,k-1)
     5 +phix(i-1,j-1,k  )
     6 +phix(i  ,j-1,k  )
     7 +phix(i-1,j  ,k  )
     8 +phix(i  ,j  ,k  ))*0.125d0
         rny00=(
     1 +phiy(i-1,j-1,k-1)
     2 +phiy(i  ,j-1,k-1)
     3 +phiy(i-1,j  ,k-1)
     4 +phiy(i  ,j  ,k-1)
     5 +phiy(i-1,j-1,k  )
     6 +phiy(i  ,j-1,k  )
     7 +phiy(i-1,j  ,k  )
     8 +phiy(i  ,j  ,k  ))*0.125d0
         rnz00=(
     1 +phiz(i-1,j-1,k-1)
     2 +phiz(i  ,j-1,k-1)
     3 +phiz(i-1,j  ,k-1)
     4 +phiz(i  ,j  ,k-1)
     5 +phiz(i-1,j-1,k  )
     6 +phiz(i  ,j-1,k  )
     7 +phiz(i-1,j  ,k  )
     8 +phiz(i  ,j  ,k  ))*0.125d0
      q00=sqrt(rnx00**2+rny00**2+rnz00**2)

      if(q00.le.eps_mthinc*dcellinv)then
      if(      phin(i,j,k).le.eps_trunc)phin(i,j,k)=0.0d0
      if(1.0d0-phin(i,j,k).le.eps_trunc)phin(i,j,k)=1.0d0
      endif
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,phin)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      phin(i,j,k)=min(phin(i,j,k),1.0d0)
      phin(i,j,k)=max(phin(i,j,k),0.0d0)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO


      sum1=0.0d0
      sum2=0.0d0
      den1=0.0d0
      den2=0.0d0
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(phi,phin)
!$OMP$ REDUCTION(+:sum1)
!$OMP$ REDUCTION(+:sum2)
!$OMP$ REDUCTION(+:den1)
!$OMP$ REDUCTION(+:den2)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      sum1=sum1+ phi(i,j,k)
      sum2=sum2+phin(i,j,k)
      den1=den1+1.0d0
      den2=den2+phin(i,j,k)**2*(1.0d0-phin(i,j,k))
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      sendbuf(1) = sum1
      sendbuf(2) = sum2
      sendbuf(3) = den1
      sendbuf(4) = den2
      recvbuf(1) = sendbuf(1)
      recvbuf(2) = sendbuf(2)
      recvbuf(3) = sendbuf(3)
      recvbuf(4) = sendbuf(4)
      call mpi_allreduce(sendbuf,recvbuf,4,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      sum1 = recvbuf(1)
      sum2 = recvbuf(2)
      den1 = recvbuf(3)
      den2 = recvbuf(4)
      endif

      offset=(sum1-sum2)/(den2+1.0d-99)
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,phin,offset)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      phin(i,j,k)=phin(i,j,k)
     & +offset*phin(i,j,k)**2*(1.0d0-phin(i,j,k))
      phin(i,j,k)=min(phin(i,j,k),1.0d0)
      phin(i,j,k)=max(phin(i,j,k),0.0d0)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

