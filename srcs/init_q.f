      subroutine init_q(ni,nj,nk,q1,q2,q3)

      implicit none
      integer ni,nj,nk,nbub
      real*8   q1(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   q2(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   q3(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(q1,q2,q3)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
        q1(i,j,k)=0.0d0
        q2(i,j,k)=0.0d0
        q3(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

        return
        end    
