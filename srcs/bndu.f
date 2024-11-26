ccc
ccc<impose boundary conditions on the velocity components uk,vk,wk
ccc<localised naier slip
ccc<l1, l2: slip length of fluid[1,2]
ccc
      subroutine bndu(nID,ni,nj,nk,uk,vk,wk,uwall,dy,l1,l2,phi)

      implicit none
      include 'param.h'
      integer nID(6)
      integer ni,nj,nk
      real*8    uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    uwall,dy,l1,l2
      real*8    phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    phi_av, ls, coef1, coef2

      integer i,j,k

ccc
ccc<j
ccc
ccc

      if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ PRIVATE(phi_av,ls,coef1,coef2)
!$OMP$ SHARED(ni,nk)
!$OMP$ SHARED(uk,vk,wk,phi,uwall,dy,l1,l2)
      do k=-2,nk+3
      do i=-2,ni+3
      phi_av = (phi(i,1,k) + phi(mod(i+1,ni+3),1,k))/2.0d0
      ls = l2 + (l2-l1)*phi_av

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      uk(i,   0,k) = coef1 * (-uwall) + coef2 * uk(i,   1,k)
      uk(i,  -1,k) = 2.d0 * uk(i,   0,k) - uk(i,   1,k)
      uk(i,  -2,k) = 2.d0 * uk(i,  -1,k) - uk(i,   0,k)

      vk(i,  -2,k)=-vk(i,   2,k)
      vk(i,  -1,k)=-vk(i,   1,k)
      vk(i,   0,k)=0.0d0

      wk(i,  -2,k)=-wk(i,   3,k)
      wk(i,  -1,k)=-wk(i,   2,k)
      wk(i,   0,k)=-wk(i,   1,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ PRIVATE(phi_av,ls,coef1,coef2)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(uk,vk,wk,phi,uwall,dy,l1,l2)
      do k=-2,nk+3
      do i=-2,ni+3
      phi_av = (phi(i,nj,k) + phi(mod(i+1,ni+3),nj,k))/2.0d0
      ls = l2 + (l2-l1)*phi_av

      coef1 = (2.d0 * dy) / (2.d0 * ls + dy)
      coef2 = (2.d0 * ls - dy) / (2.d0 * ls + dy)

      uk(i,nj+1,k) = coef1 * uwall + coef2 * uk(i,nj  ,k)
      uk(i,nj+2,k) = 2.d0 * uk(i,nj+1,k) - uk(i,nj  ,k)
      uk(i,nj+3,k) = 2.d0 * uk(i,nj+2,k) - uk(i,nj+1,k)

      vk(i,nj  ,k)=0.0d0
      vk(i,nj+1,k)=-vk(i,nj-1,k)
      vk(i,nj+2,k)=-vk(i,nj-2,k)
      vk(i,nj+3,k)=-vk(i,nj-3,k)

      wk(i,nj+1,k)=-wk(i,nj  ,k)
      wk(i,nj+2,k)=-wk(i,nj-1,k)
      wk(i,nj+3,k)=-wk(i,nj-2,k)
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      return
      end

