      subroutine caldt(ipara,nID,ID,ndiv,ni,nj,nk,nstep,imon_t
     & ,dxinv,dyinv,dzinv
     & ,cfl,rhol,rhog,rmul,rmug,surface_tension
     & ,u,v,w
     & ,dt,time)

      implicit none
      include 'mpif.h'
      include 'param.h'

      integer ipara,nID(6),ID,ndiv
      integer ni,nj,nk
      integer nstep,imon_t
      real*8 dxinv,dyinv,dzinv
      real*8 cfl,rhol,rhog,rmul,rmug,surface_tension
      real*8 u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 dt,time

      real*8 recvbuf(1),sendbuf(1)
      integer ierr
      integer jen
      integer i,j,k
      real*8 advmax,dcell,rnu,st_rhominv,dt1,dt2,dt3

      dcell=(dxinv*dyinv*dzinv)**(-1.0d0/3.0d0)
      rnu=max(rmul/rhol,rmug/rhog)
      st_rhominv=2.0d0*surface_tension/(rhol+rhog)

      jen=nj
      if(nID(Y_PLUS).lt.0)jen=nj-1

ccc
ccc<computate maximum advection speed
ccc
      advmax=1.0d-20

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(u,dxinv)
!$OMP$ REDUCTION(max:advmax)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      advmax=max(abs(u(i,j,k)*dxinv),advmax)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,jen,nk)
!$OMP$ SHARED(v,dyinv)
!$OMP$ REDUCTION(max:advmax)
      do k=1,nk
      do j=1,jen
      do i=1,ni
      advmax=max(abs(v(i,j,k)*dyinv),advmax)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(w,dzinv)
!$OMP$ REDUCTION(max:advmax)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      advmax=max(abs(w(i,j,k)*dzinv),advmax)
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      if(ipara.eq.1)then
      sendbuf(1) = advmax
      recvbuf(1) = advmax
      call mpi_allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_MAX
     & ,MPI_COMM_WORLD,ierr)
      advmax = recvbuf(1)
      endif

ccc
ccc<detemine time increament dt and update time
ccc dt1: due to advection
ccc dt2: due to viscosity
ccc dt3: due to surface tension
ccc

      dt1=cfl/advmax
      dt2=cfl*0.5d0*dcell**2/rnu*1.0d2
      dt3=cfl*sqrt(dcell**3/st_rhominv)

      dt=min(dt1,dt2)
      dt=min(dt ,dt3)

      time=time+dt
      if(mod(nstep,imon_t).eq.0.and.ID.eq.0)then
      write(*,'("time=",1e17.10," dt=",1e17.10
     &         ," dt_adv=",1e17.10
     &         ," dt_vis=",1e17.10
     &         ," dt_srf=",1e17.10)')
     & time,dt,dt1,dt2,dt3
      endif

      return
      end

