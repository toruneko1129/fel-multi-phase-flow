      subroutine cal_grad_p2uvw(ID,svall,ni,nj,nk
     & ,dxinv,dyinv_array,dzinv
     & ,q,qx_u,qy_v,qz_w)

      implicit none
      integer ni,nj,nk,ID,svall
      real*8 dxinv,dyinv,dzinv
      real*8 dyinv_array(-2:svall+3)
      real*8    q(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qx_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qy_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 qz_w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k,dis
      real*8 q_000,q_p00,q_0p0,q_00p

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,dis,dyinv)
!$OMP$ PRIVATE(q_000,q_p00,q_0p0,q_00p)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(q,qx_u,qy_v,qz_w)
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
      q_00p=q(i  ,j  ,k+1)

      qx_u(i,j,k)=(-q_000+q_p00)*dxinv
      qy_v(i,j,k)=(-q_000+q_0p0)*dyinv
      qz_w(i,j,k)=(-q_000+q_00p)*dzinv
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

