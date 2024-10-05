      subroutine cal_vel(ipara,ni,nj,nk,radius,rhol,rmul
     & ,u,v,w,phi,vel,re)

      implicit none
      integer nID(6)
      include 'mpif.h'
      integer ipara,ni,nj,nk
      real*8     radius,rhol,rmul
      real*8     phi(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8     w(-2:ni+3,-2:nj+3,-2:nk+3)

      integer ierr
      real*8  sendbuf(4),recvbuf(4)
      integer i,j,k
      real*8  delta,sum_phi,vel(3),re(3)
      real*8  phi_u,phi_v,phi_w,diameter
      
      vel(1)=0.0d0
      vel(2)=0.0d0
      vel(3)=0.0d0
      re(1) =0.0d0
      re(2) =0.0d0
      re(3) =0.0d0
      sum_phi=0.0d0
      phi_u=0.0d0     
      phi_v=0.0d0
      phi_w=0.0d0
      diameter=radius*2.0d0
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(phi,u,v,w)
!$OMP$ REDUCTION(+:sum_phi)
!$OMP$ REDUCTION(+:phi_u)
!$OMP$ REDUCTION(+:phi_v)
!$OMP$ REDUCTION(+:phi_w)

      !calculate the summation of phi (donominator)
      !calculate the summation of phi*u (numerator)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sum_phi=sum_phi+phi(i,j,k)
        
        phi_u=phi_u+phi(i,j,k)*u(i,j,k)
        phi_v=phi_v+phi(i,j,k)*v(i,j,k)
        phi_w=phi_w+phi(i,j,k)*w(i,j,k)
      enddo
      enddo
      enddo
      
!$OMP  END PARALLEL DO
      
      if(ipara.eq.1)then
      sendbuf(1) = sum_phi
      sendbuf(2) = phi_u
      sendbuf(3) = phi_v
      sendbuf(4) = phi_w
      recvbuf(1) = sendbuf(1)
      recvbuf(2) = sendbuf(2)
      recvbuf(3) = sendbuf(3)
      recvbuf(4) = sendbuf(4)
      call mpi_allreduce(sendbuf,recvbuf,4,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      sum_phi= recvbuf(1)
      phi_u  = recvbuf(2)
      phi_v  = recvbuf(3)
      phi_w  = recvbuf(4)
      endif

      !calculate the velocity of bubble
      vel(1)=phi_u/sum_phi
      vel(2)=phi_v/sum_phi
      vel(3)=phi_w/sum_phi
      re(1) =rhol*vel(1)*diameter/rmul
      re(2) =rhol*vel(2)*diameter/rmul
      re(3) =rhol*vel(3)*diameter/rmul
      !write(*,'("vel",3e20.13)')vel(1),vel(2),vel(3)
      return
      end
