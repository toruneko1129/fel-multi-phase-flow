ccc
ccc solve pressure equation using FFT and TDMA
ccc
      subroutine solp_fft_tdma3(ni,nj,nk
     & ,fftdata,phir_r,phir_i,dp)

      implicit none
      integer ni,nj,nk
      real*8  fftdata(ni*nk*2,nj)
      real*8   phir_r(-ni/2:ni/2-1,nj,-nk/2:nk/2-1)
      real*8   phir_i(-ni/2:ni/2-1,nj,-nk/2:nk/2-1)
      real*8       dp(-2:ni+3,-2:nj+3,-2:nk+3)

      integer nn(2)
      integer i,j,k,l,ii,kk

      nn(1)=ni
      nn(2)=nk

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,ii,kk,l)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(fftdata,phir_r,phir_i,nn,dp)
      do j=1,nj
      do k=-nk/2,nk/2-1
      kk=mod(k+nk,nk)
      do i=-ni/2,ni/2-1
      ii=mod(i+ni,ni)
      l=ni*kk+ii+1
      fftdata(2*l-1,j)=phir_r(i,j,k)
      fftdata(2*l  ,j)=phir_i(i,j,k)
      enddo
      enddo

      call fourn(fftdata(1,j),nn,2,1)

      do k=0,nk-1
      do i=0,ni-1
      l=ni*k+i+1
      dp(i+1,j,k+1)=fftdata(2*l-1,j)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

c -------------------------------------------------------------

