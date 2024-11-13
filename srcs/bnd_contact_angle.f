ccc
ccc<impose the contact angle boundary condition on qk
ccc<give the contact angle theta_deg and grid space dy
ccc
      subroutine bnd_contact_angle(nID, ni, nj, nk, qk, theta_deg, dy)

      implicit none
      include 'param.h'
      integer nID(6),ni,nj,nk
      real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 theta_deg, dy

      integer i,k
      real*8 theta_rad, cot_theta


      ! Convert angle from degrees to radians
      theta_rad = theta_deg * (acos(-1.0d0) / 180.0d0)
      cot_theta = 1.0d0 / tan(theta_rad)

      if(nID(Y_MINUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nk,qk,cot_theta,dy)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i, -2, k) = qk(i, 3, k) + 5 * dy * cot_theta
      qk(i, -1, k) = qk(i, 2, k) + 3 * dy * cot_theta
      qk(i,  0, k) = qk(i, 1, k) +     dy * cot_theta
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      if(nID(Y_PLUS).lt.0)then
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(ni,nj,nk,qk,cot_theta,dy)
      do k=-2,nk+3
      do i=-2,ni+3
      qk(i, nj + 1, k) = qk(i, nj    , k) +     dy * cot_theta
      qk(i, nj + 2, k) = qk(i, nj - 1, k) + 3 * dy * cot_theta
      qk(i, nj + 3, k) = qk(i, nj - 2, k) + 5 * dy * cot_theta
      enddo
      enddo
!$OMP  END PARALLEL DO
      endif

      return
      end
