      subroutine cal_grad_p2a(ID,svall,ni,nj,nk
     & ,dxinv,dyinv_array,dzinv,q,qx,qy,qz)

      implicit none
      integer ni,nj,nk,ID,svall
      real*8 dxinv,dyinv,dzinv
      real*8  q(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qx(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qy(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qz(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 dyinv_array(-2:svall+3)

      integer i,j,k,dis
      real*8 q_000,q_p00,q_0p0,q_pp0
      real*8 q_00p,q_p0p,q_0pp,q_ppp
      real*8 norm,norminv,eps

      eps = 1.0d-20

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,dis,dyinv)
!$OMP$ PRIVATE(q_000,q_p00,q_0p0,q_pp0)
!$OMP$ PRIVATE(q_00p,q_p0p,q_0pp,q_ppp)
!$OMP$ PRIVATE(norm,norminv)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(q,qx,qy,qz,eps)
!$OMP$ SHARED(dxinv,dyinv_array,dzinv)
      do k=-2,nk+2
      do j=-2,nj+2
      do i=-2,ni+2
      
      dis=j+nj*ID
      !mesh number from the wall
      dyinv=dyinv_array(dis)
      q_000=q(i  ,j  ,k  )
      q_p00=q(i+1,j  ,k  )
      q_0p0=q(i  ,j+1,k  )
      q_pp0=q(i+1,j+1,k  )
      q_00p=q(i  ,j  ,k+1)
      q_p0p=q(i+1,j  ,k+1)
      q_0pp=q(i  ,j+1,k+1)
      q_ppp=q(i+1,j+1,k+1)

      qx(i,j,k)=-q_000+q_p00-q_0p0+q_pp0-q_00p+q_p0p-q_0pp+q_ppp
      qy(i,j,k)=-q_000-q_p00+q_0p0+q_pp0-q_00p-q_p0p+q_0pp+q_ppp
      qz(i,j,k)=-q_000-q_p00-q_0p0-q_pp0+q_00p+q_p0p+q_0pp+q_ppp

      !normalise
      !norm = sqrt(qx(i,j,k)**2 + qy(i,j,k)**2 + qz(i,j,k)**2)
      !norminv = 1.0d0/max(norm,eps);
      
      !qx(i,j,k)=qx(i,j,k)*norminv
      !qy(i,j,k)=qy(i,j,k)*norminv
      !qz(i,j,k)=qz(i,j,k)*norminv

      qx(i,j,k)=qx(i,j,k)*0.25d0*dxinv
      qy(i,j,k)=qy(i,j,k)*0.25d0*dyinv
      qz(i,j,k)=qz(i,j,k)*0.25d0*dzinv
      
      !to define just the location of wall
      !if(dis.eq.0)then
      !qx(i,j,k)=qx(i,j,k)*2.0d0
      !qy(i,j,k)=qy(i,j,k)*2.0d0
      !qz(i,j,k)=qz(i,j,k)*2.0d0
      !endif

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      !if(norm.lt.eps)write(*,'("wrong")')
      return
      end

