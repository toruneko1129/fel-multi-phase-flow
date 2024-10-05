      subroutine cal_center(ipara,ID,svall,ni,nj,nk,xl,yl,zl,dt
     & ,rhog,phi,center,center_pre1,center_pre2
     & ,area_pre1,area_pre2,lap,velocity)
      implicit none
      integer nID(6)
      include 'mpif.h'
      integer ipara,ID,svall,ni,nj,nk
      real*8     xl,yl,zl,rhog,dt
      real*8     phi(-2:ni+3,-2:nj+3,-2:nk+3)
      integer area_pre1,area_pre2,lap
      real*8     center_pre1,center_pre2

      integer ierr
      real*8 sendbuf(4),recvbuf(4)
      integer i,j,k
      real*8  delta,weight,volume,center(3),center_cal,location(3)
      real*8  we_lo_x,we_lo_y,we_lo_z,sensor1,sensor2,velocity
      
      center(1)=0.0d0
      center(2)=0.0d0
      center(3)=0.0d0
      weight=0.0d0
      we_lo_x=0.0d0
      we_lo_y=0.0d0
      we_lo_z=0.0d0
      delta=yl/svall
      volume=delta*delta*delta  !volume: per 1 mesh
      sensor1=xl*0.25d0
      sensor2=xl*0.75d0
      !detect how many laps bubbles move
      if(area_pre1 .eq. 1 .and. area_pre2 .eq. 3)lap=lap+1

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(location)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(phi,xl,delta,lap,area_pre1,area_pre2,sensor1,sensor2)
!$OMP$ REDUCTION(+:weight)
!$OMP$ REDUCTION(+:we_lo_x)
!$OMP$ REDUCTION(+:we_lo_y)
!$OMP$ REDUCTION(+:we_lo_z)

      !calculate the total weight:weight (donominator)
      !calculate the weight*location:we_lo (numerator)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        weight=weight+phi(i,j,k)
        
        location(1)=delta*(i-0.5d0)
        location(2)=delta*(j+ID*nj-0.5d0)
        location(3)=delta*(k-0.5d0)

        !Area 1,2,3  bubble location
        !1: 0-sensor1
        !2: sensor1-sensor2
        !3: sensor2-xl
        if(area_pre1.eq.3 .and. location(1).lt.sensor1)then
          location(1)=location(1)+xl
        endif
        if(area_pre1.eq.1 .and. location(1).gt.sensor2)then
          location(1)=location(1)-xl
        endif

        location(1)=location(1)+xl*lap
        we_lo_x=we_lo_x+phi(i,j,k)*location(1)
        we_lo_y=we_lo_y+phi(i,j,k)*location(2)
        we_lo_z=we_lo_z+phi(i,j,k)*location(3)
      enddo
      enddo
      enddo
      
!$OMP  END PARALLEL DO
      
      if(ipara.eq.1)then
      sendbuf(1) = weight
      sendbuf(2) = we_lo_x
      sendbuf(3) = we_lo_y
      sendbuf(4) = we_lo_z
      recvbuf(1) = sendbuf(1)
      recvbuf(2) = sendbuf(2)
      recvbuf(3) = sendbuf(3)
      recvbuf(4) = sendbuf(4)
      call mpi_allreduce(sendbuf,recvbuf,4,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      weight = recvbuf(1)
      we_lo_x= recvbuf(2)
      we_lo_y= recvbuf(3)
      we_lo_z= recvbuf(4)
      endif

      !calculate the center of gravity of bubble
      center(1)=we_lo_x/weight
      center(2)=we_lo_y/weight
      center(3)=we_lo_z/weight
      !write(*,'("center1",3e20.13)')center(1),center(2),center(3)

      !calculate the velocity
      velocity    = (center(1)-center_pre1)/dt

      center_pre2 = center_pre1
      center_pre1 = center(1)
      center_cal  = center(1)
      !detect the area(1,2,3) of bubbles in previous step and this step
      area_pre2   = area_pre1
      !Repeat
      if(center_cal.gt.xl)then
      do
      center_cal=center_cal-xl
      if(center_cal.lt.xl)exit
      enddo
      endif

      if(center_cal.lt.sensor1)area_pre1=1
      if(center_cal.gt.sensor2)area_pre1=3
      if(sensor1.le.center_cal .and. center_cal.le.sensor2)area_pre1=2

      return
      end
