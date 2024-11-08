subroutine init_phi(ndiv, svall, phi, nbub, num, nsv)

      implicit none
      integer ndiv, svall(3), nbub, num
      integer nsv
      real*8 phi(-2:svall(1)+3, -2:svall(2)/ndiv+3, -2:svall(3)+3, 0:nbub)

      integer i, j, k
      integer ni, nj, nk

      ni = svall(1)
      nj = svall(2) / ndiv
      nk = svall(3)

!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i, j, k) SHARED(ni, nj, nk, phi, num, nsv)
      do k = -2, nk + 3
      do j = 1, nj
      do i = -2, ni + 3
          if (2 * nsv < i .and. i <= 6 * nsv) then
              phi(i, j, k, num) = 1.0d0
          else
              phi(i, j, k, num) = 0.0d0
          endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      return
      end