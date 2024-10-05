ccc
ccc solve pressure equation using FFT and TDMA
ccc
      subroutine solp_fft_tdma2(ID,ndiv,svall
     & ,rho,dxinv,dzinv
     & ,as_p_fft,an_p_fft,ap_p_fft,atdma,btdma_r,btdma_i
     & ,phiw_r,phiw_i)

      implicit none
      integer ID,ndiv,svall(3)
      real*8 rho,dxinv,dzinv
      real*8 as_p_fft(-2:svall(2)+3)
      real*8 an_p_fft(-2:svall(2)+3)
      real*8 ap_p_fft(-2:svall(2)+3)
      real*8    atdma(-2:svall(2)+3)
      real*8  btdma_r(-2:svall(2)+3)
      real*8  btdma_i(-2:svall(2)+3)
      real*8   phiw_r(-svall(1)/2:svall(1)/2-1,svall(2),svall(3)/ndiv)
      real*8   phiw_i(-svall(1)/2:svall(1)/2-1,svall(2),svall(3)/ndiv)

      integer ni,nj,nk
      integer i,j,k
      real*8 pi,dxinv2_rho,dzinv2_rho,pi2_ni,pi2_nk
      real*8 thtx,thtz,c1x,c1z,ap,d00,apinv

      ni=svall(1)
      nj=svall(2)
      nk=svall(3)/ndiv

      pi=atan(1.0d0)*4.0d0
      dxinv2_rho=dxinv**2/rho
      dzinv2_rho=dzinv**2/rho
      pi2_ni =pi*2.0d0/dble(svall(1))
      pi2_nk =pi*2.0d0/dble(svall(3))

!$OMP PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(thtx,c1x,thtz,c1z)
!$OMP$ PRIVATE(atdma,btdma_r,btdma_i)
!$OMP$ PRIVATE(ap,d00,apinv)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(ID,svall,pi2_ni,pi2_nk,dxinv2_rho,dzinv2_rho)
!$OMP$ SHARED(ap_p_fft,as_p_fft,an_p_fft)
!$OMP$ SHARED(phiw_r,phiw_i)
      do k=1,nk
      do i=-ni/2,ni/2-1
      thtx=dble(i)*pi2_ni
      c1x=cos(thtx)*2.0d0
      thtz=dble(k+ID*nk-svall(3)/2-1)*pi2_nk
      c1z=cos(thtz)*2.0d0

c
c<step 1
c<neumann condition in y direction
        atdma(   0)=0.0d0
      btdma_r(   0)=0.0d0
      btdma_i(   0)=0.0d0
        atdma(nj+1)=0.0d0
      btdma_r(nj+1)=0.0d0
      btdma_i(nj+1)=0.0d0

c
c<step 2
c
      do j=1,nj
      ap=ap_p_fft(j)
     &   +dxinv2_rho*(2.0d0-c1x)
     &   +dzinv2_rho*(2.0d0-c1z)
      d00=-ap +as_p_fft(j)*atdma(j-1)
      apinv=-d00/(d00**2+1.0d-99)
        atdma(j)=apinv*an_p_fft(j)
      btdma_r(j)=apinv*(-phiw_r(i,j,k)+as_p_fft(j)*btdma_r(j-1))
      btdma_i(j)=apinv*(-phiw_i(i,j,k)+as_p_fft(j)*btdma_i(j-1))
      enddo

c
c<step 3
c
      phiw_r(i,nj,k)=btdma_r(nj)
      phiw_i(i,nj,k)=btdma_i(nj)
      do j=nj-1,1,-1
      phiw_r(i,j,k)=btdma_r(j)+atdma(j)*phiw_r(i,j+1,k)
      phiw_i(i,j,k)=btdma_i(j)+atdma(j)*phiw_i(i,j+1,k)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

c -------------------------------------------------------------

