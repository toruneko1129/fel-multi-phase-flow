      subroutine summation(ni,nj,nk,qk,nbub)

      implicit none
      integer ni,nj,nk,nbub
      real*8 qk(-2:ni+3,-2:nj+3,-2:nk+3,0:nbub)

      integer i,j,k,l

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,qk)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      qk(i,j,k,0)=qk(i,j,k,1)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

      do l=2,nbub
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(l,ni,nj,nk,nbub,qk)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      qk(i,j,k,0)=qk(i,j,k,0)+qk(i,j,k,l)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      enddo

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,nbub,qk)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      qk(i,j,k,0)=min(max(qk(i,j,k,0),0.0d0),1.0d0)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      
      return
      end
