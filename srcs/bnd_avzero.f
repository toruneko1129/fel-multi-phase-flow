      subroutine bnd_avzero(ipara,ndiv,ni,nj,nk,q)

      implicit none
      integer ipara,ndiv,ni,nj,nk
      real*8 q(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k
      real*8 av

      call cal_av(ipara,ndiv,ni,nj,nk,q,av)

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(q,av)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      q(i,j,k)=q(i,j,k)-av
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end

