ccc
ccc solve pressure equation using FFT and TDMA
ccc
      subroutine solp_fft_tdma1(ni,nj,nk
     & ,fftdata,div,phir_r,phir_i)

      implicit none
      integer ni,nj,nk
      real*8  fftdata(ni*nk*2,nj)
      real*8      div(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   phir_r(-ni/2:ni/2-1,nj,-nk/2:nk/2-1)
      real*8   phir_i(-ni/2:ni/2-1,nj,-nk/2:nk/2-1)

      integer nn(2)
      integer i,j,k,l,ii,kk

      nn(1)=ni
      nn(2)=nk

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,l,ii,kk)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(fftdata,div,nn,phir_r,phir_i)
      do j=1,nj

      do k=0,nk-1
      do i=0,ni-1
      l=ni*k+i+1
      fftdata(2*l-1,j)=div(i+1,j,k+1)
      fftdata(2*l  ,j)=0.0d0
      enddo
      enddo

      call fourn(fftdata(1,j),nn,2,-1)

      do k=-nk/2,nk/2-1
      kk=mod(k+nk,nk)
      do i=-ni/2,ni/2-1
      ii=mod(i+ni,ni)
      l=ni*kk+ii+1
      phir_r(i,j,k)=fftdata(2*l-1,j)/dble(ni*nk)
      phir_i(i,j,k)=fftdata(2*l  ,j)/dble(ni*nk)
      enddo
      enddo

      enddo
!$OMP  END PARALLEL DO
      return
      end
