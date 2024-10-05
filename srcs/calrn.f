      subroutine calrn(ni,nj,nk,dxinv,dyinv,dzinv
     & ,ip,jp,kp,phix,phiy,phiz,rn)

      implicit none
      integer ni,nj,nk
      integer ip,jp,kp
      real*8 dxinv,dyinv,dzinv
      real*8  phix(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  phiy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  phiz(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    rn(-2:ni+3,-2:nj+3,-2:nk+3,3)

      integer i,j,k,ii,jj,kk
      real*8 rnx00,rny00,rnz00,arn00,arninv,verysmall

      verysmall=sqrt(dxinv**2+dyinv**2+dzinv**2)*1.0d-30

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,rn)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      rn(i,j,k,1)=0.0d0
      rn(i,j,k,2)=0.0d0
      rn(i,j,k,3)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO


ccc
ccc<compute normal vector rn
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,ii,jj,kk)
!$OMP$ PRIVATE(rnx00,rny00,rnz00,arn00,arninv)
!$OMP$ SHARED(ni,nj,nk,ip,jp,kp,phix,phiy,phiz,verysmall,rn)
      do k=-2,nk+3-kp
      do j=-2,nj+3-jp
      do i=-2,ni+3-ip
      ii=i+ip
      jj=j+jp
      kk=k+kp
         rnx00=(
     1 +phix(i ,j ,k )
     2 +phix(ii,j ,k )
     3 +phix(i ,jj,k )
     4 +phix(ii,jj,k )
     5 +phix(i ,j ,kk)
     6 +phix(ii,j ,kk)
     7 +phix(i ,jj,kk)
     8 +phix(ii,jj,kk))*0.125d0
         rny00=(
     1 +phiy(i ,j ,k )
     2 +phiy(ii,j ,k )
     3 +phiy(i ,jj,k )
     4 +phiy(ii,jj,k )
     5 +phiy(i ,j ,kk)
     6 +phiy(ii,j ,kk)
     7 +phiy(i ,jj,kk)
     8 +phiy(ii,jj,kk))*0.125d0
         rnz00=(
     1 +phiz(i ,j ,k )
     2 +phiz(ii,j ,k )
     3 +phiz(i ,jj,k )
     4 +phiz(ii,jj,k )
     5 +phiz(i ,j ,kk)
     6 +phiz(ii,j ,kk)
     7 +phiz(i ,jj,kk)
     8 +phiz(ii,jj,kk))*0.125d0

      arn00=sqrt(rnx00**2+rny00**2+rnz00**2)+verysmall
      arninv=1.0d0/arn00

      rn(ii,jj,kk,1)=rnx00*arninv
      rn(ii,jj,kk,2)=rny00*arninv
      rn(ii,jj,kk,3)=rnz00*arninv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end




