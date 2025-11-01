ccc
ccc<impose boundary conditions on the velocity components uk,vk,wk
ccc<localised naier slip
ccc<l1, l2: slip length of fluid[1,2]
ccc
      subroutine bndu(nID,ni,nj,nk,uk,vk,wk,u_top,u_bot,dy,l1,l2,phi)

      implicit none
      include 'param.h'
      integer nID(6)
      integer ni,nj,nk
      real*8    uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    u_top,u_bot,dy
      real*8    l1(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    l2(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    inv_ls1, inv_ls2
      real*8    phi_av, ls, coef1, coef2, eps

      integer i,j,k,ip,kp

ccc
ccc<j
ccc
ccc

      eps=1.0d-12

      if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k,ip,kp)
!$OMP$ PRIVATE(phi_av,inv_ls1,inv_ls2,ls,coef1,coef2)
!$OMP$ SHARED(ni,nk)
!$OMP$ SHARED(uk,vk,wk,phi,u_bot,dy,l1,l2,eps)
      do k=1,nk
      do i=1,ni

      ip = i + 1
      if (ip .gt. ni) ip = 1
      kp = k + 1
      if (kp .gt. nk) kp = 1

      phi_av = (phi(i,1,k) + phi(ip,1,k))/2.0d0
      inv_ls1 = phi_av / (l1(i,1,k)+eps)
      inv_ls2 = (1.d0 - phi_av) / (l2(i,1,k)+eps)
      ls = 1.d0 / (inv_ls1 + inv_ls2+eps)

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      uk(i,   0,k) = coef1 * u_bot + coef2 * uk(i,   1,k)
      uk(i,  -1,k) = 2.d0 * uk(i,   0,k) - uk(i,   1,k)
      uk(i,  -2,k) = 2.d0 * uk(i,  -1,k) - uk(i,   0,k)

      vk(i,  -2,k)=-vk(i,   2,k)
      vk(i,  -1,k)=-vk(i,   1,k)
      vk(i,   0,k)=0.0d0

      phi_av = (phi(i,1,k) + phi(i,1,kp))/2.0d0
      inv_ls1 = phi_av / (l1(i,1,k)+eps)
      inv_ls2 = (1.d0 - phi_av) / (l2(i,1,k)+eps)
      ls = 1.d0 / (inv_ls1 + inv_ls2+eps)

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      wk(i,   0,k) = coef2 * wk(i,   1,k)
      wk(i,  -1,k) = 2.d0 * wk(i,   0,k) - wk(i,   1,k)
      wk(i,  -2,k) = 2.d0 * wk(i,  -1,k) - wk(i,   0,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k,ip,kp)
!$OMP$ PRIVATE(phi_av,inv_ls1,inv_ls2,ls,coef1,coef2)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,vk,wk,phi,u_top,dy,l1,l2,eps)
      do k=1,nk
      do i=1,ni

      ip = i + 1
      if (ip .gt. ni) ip = 1
      kp = k + 1
      if (kp .gt. nk) kp = 1

      phi_av = (phi(i,nj,k) + phi(ip,nj,k))/2.0d0
      inv_ls1 = phi_av / (l1(i,nj,k)+eps)
      inv_ls2 = (1.d0 - phi_av) / (l2(i,nj,k)+eps)
      ls = 1.d0 / (inv_ls1 + inv_ls2+eps)

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      uk(i,nj+1,k) = coef1 * u_top + coef2 * uk(i,nj  ,k)
      uk(i,nj+2,k) = 2.d0 * uk(i,nj+1,k) - uk(i,nj  ,k)
      uk(i,nj+3,k) = 2.d0 * uk(i,nj+2,k) - uk(i,nj+1,k)

      vk(i,nj  ,k)=0.0d0
      vk(i,nj+1,k)=-vk(i,nj-1,k)
      vk(i,nj+2,k)=-vk(i,nj-2,k)
      vk(i,nj+3,k)=-vk(i,nj-3,k)

      phi_av = (phi(i,nj,k) + phi(i,nj,kp))/2.0d0
      inv_ls1 = phi_av / (l1(i,nj,k)+eps)
      inv_ls2 = (1.d0 - phi_av) / (l2(i,nj,k)+eps)
      ls = 1.d0 / (inv_ls1 + inv_ls2+eps)

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      wk(i,nj+1,k) = coef2 * wk(i,nj  ,k)
      wk(i,nj+2,k) = 2.d0 * wk(i,nj+1,k) - wk(i,nj  ,k)
      wk(i,nj+3,k) = 2.d0 * wk(i,nj+2,k) - wk(i,nj+1,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,vk,wk)

      do j = 1, nj
      do k = 1, nk
      uk(0 ,j,k) = uk(ni+1,j,k)
      uk(-1,j,k) = uk(ni  ,j,k)
      uk(-2,j,k) = uk(ni-1,j,k)
      vk(0 ,j,k) = vk(ni  ,j,k)
      vk(-1,j,k) = vk(ni-1,j,k)
      vk(-2,j,k) = vk(ni-2,j,k)
      wk(0 ,j,k) = wk(ni  ,j,k)
      wk(-1,j,k) = wk(ni-1,j,k)
      wk(-2,j,k) = wk(ni-2,j,k)

      uk(ni+2, j, k) = uk(2, j, k)
      uk(ni+3, j, k) = uk(3, j, k)
      vk(ni+1, j, k) = vk(1, j, k)
      vk(ni+2, j, k) = vk(2, j, k)
      vk(ni+3, j, k) = vk(3, j, k)
      wk(ni+1, j, k) = wk(1, j, k)
      wk(ni+2, j, k) = wk(2, j, k)
      wk(ni+3, j, k) = wk(3, j, k)
      end do
      end do

!$OMP  END PARALLEL DO


!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j)
!$OMP$ SHARED(ni,nj,nk,uk,vk,wk)
      do i = 1, ni
      do j = 1, nj
      uk(i, j,   0) = uk(i, j, nk  )
      uk(i, j,  -1) = uk(i, j, nk-1)
      uk(i, j,  -2) = uk(i, j, nk-2)
      vk(i, j,   0) = vk(i, j, nk  )
      vk(i, j,  -1) = vk(i, j, nk-1)
      vk(i, j,  -2) = vk(i, j, nk-2)
      wk(i, j,   0) = wk(i, j, nk+1)
      wk(i, j,  -1) = wk(i, j, nk  )
      wk(i, j,  -2) = wk(i, j, nk-1)

      uk(i, j, nk+1) = uk(i, j,   1)
      uk(i, j, nk+2) = uk(i, j,   2)
      uk(i, j, nk+3) = uk(i, j,   3)
      vk(i, j, nk+1) = vk(i, j,   1)
      vk(i, j, nk+2) = vk(i, j,   2)
      vk(i, j, nk+3) = vk(i, j,   3)
      wk(i, j, nk+2) = wk(i, j,   2)
      wk(i, j, nk+3) = wk(i, j,   3)
      end do
      end do
!$OMP  END PARALLEL DO

      return
      end

