ccc
ccc<compute the viscous stress components: tau[ij]k
ccc from the strain rate components s[ij]k
ccc and the viscosity rmuk
ccc

      subroutine cal_arith_tau(ni,nj,nk,s,rmu,tau)
      implicit none
      integer ni,nj,nk
      real*8     s(-2:ni+3,-2:nj+3,-2:nk+3,6)
      real*8   rmu(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   tau(-2:ni+3,-2:nj+3,-2:nk+3,6)

      integer i,j,k,l
      real*8 rmu00

ccc
ccc<for diagonal components:
ccc

      do l=1,3
!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk,l)
!$OMP$ SHARED(tau,s,rmu)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
      tau(i,j,k,l)=s(i,j,k,l)*rmu(i,j,k)*2.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      enddo

ccc
ccc<for non-diagonal components:
ccc

ccc
ccc \tau_{12} = 2\mu S_{12}
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rmu00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(tau,s,rmu)
      do k=-2,nk+3
      do j=-2,nj+2
      do i=-2,ni+2

      rmu00=
     1 +0.25d0*rmu(i  ,j  ,k  )
     2 +0.25d0*rmu(i+1,j  ,k  )
     3 +0.25d0*rmu(i  ,j+1,k  )
     4 +0.25d0*rmu(i+1,j+1,k  )

      tau(i,j,k,4)=s(i,j,k,4)*rmu00*2.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc \tau_{13} = 2\mu S_{13}
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rmu00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(tau,s,rmu)
      do k=-2,nk+2
      do j=-2,nj+3
      do i=-2,ni+2

      rmu00=
     1 +0.25d0*rmu(i  ,j  ,k  )
     2 +0.25d0*rmu(i+1,j  ,k  )
     3 +0.25d0*rmu(i  ,j  ,k+1)
     4 +0.25d0*rmu(i+1,j  ,k+1)

      tau(i,j,k,5)=s(i,j,k,5)*rmu00*2.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc \tau_{23} = 2\mu S_{23}
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k,rmu00)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(tau,s,rmu)
      do k=-2,nk+2
      do j=-2,nj+2
      do i=-2,ni+3

      rmu00=
     1 +0.25d0*rmu(i  ,j  ,k  )
     2 +0.25d0*rmu(i  ,j+1,k  )
     3 +0.25d0*rmu(i  ,j  ,k+1)
     4 +0.25d0*rmu(i  ,j+1,k+1)

      tau(i,j,k,6)=s(i,j,k,6)*rmu00*2.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      
      return
      end
