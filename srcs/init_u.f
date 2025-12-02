      subroutine init_u(ni,nj,nk,q1,u_in)

      implicit none
      integer ni,nj,nk
      real*8   q1(-2:ni+3,-2:nj+3,-2:nk+3)
	  real*8   u_in

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u_in,q1)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
        q1(i,j,k)=u_in
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

        return
        end    
