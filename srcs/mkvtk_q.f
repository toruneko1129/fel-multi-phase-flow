      subroutine mkvtk_q(svall,nstep,dx,dy,dz,vorx,q)
      implicit none
      integer svall(3),nstep
      real*8 dx,dy,dz
      real*8 vorx(-2:svall(1)+3,-2:svall(2)+3,-2:svall(3)+3)
      real*8    q(-2:svall(1)+3,-2:svall(2)+3,-2:svall(3)+3)

      character*32 fname
      integer i,j,k,mi,mj,mk
      real*8 vorx00,q00

      mi=svall(1)
      mj=svall(2)
      mk=svall(3)

      write(fname,'("q1_",i7.7,".vtk")')nstep
      open(10,file=fname,status='replace')
      write(10,"('# vtk DataFile Version 2.0')")
      write(10,"('t=',I0)")nstep
      write(10,"('ASCII')")
      write(10,"('DATASET STRUCTURED_POINTS')")
      write(10,"('DIMENSIONS ',I0,' ',I0,' ',I0)")mi,mj,mk
      write(10,"('ORIGIN ',e11.4,' ',e11.4,' ',e11.4)")
     & dx*0.5d0,dy*0.5d0,dz*0.5d0
      write(10,"('ASPECT_RATIO ',e11.4,' ',e11.4,' ',e11.4)")dx,dy,dz
      write(10,"('')")
      write(10,"('POINT_DATA ',I0)")mi*mj*mk
      write(10,"('SCALARS p float')")
      write(10,"('LOOKUP_TABLE default')")
      do k=1,mk
      do j=1,mj
      do i=1,mi
      vorx00=(
     1  vorx(i,j-1,k-1)
     2 +vorx(i,j  ,k-1)
     3 +vorx(i,j-1,k  )
     4 +vorx(i,j  ,k  ))*0.25d0
      q00=q(i,j,k)
      if(vorx00.le.0.0d0)q00=0.0d0
      write(10,"(e10.3)")q00
      enddo
      enddo
      enddo
      close(10)
      write(*,*)'output ',fname

      write(fname,'("q2_",i7.7,".vtk")')nstep
      open(10,file=fname,status='replace')
      write(10,"('# vtk DataFile Version 2.0')")
      write(10,"('t=',I0)")nstep
      write(10,"('ASCII')")
      write(10,"('DATASET STRUCTURED_POINTS')")
      write(10,"('DIMENSIONS ',I0,' ',I0,' ',I0)")mi,mj,mk
      write(10,"('ORIGIN ',e11.4,' ',e11.4,' ',e11.4)")
     & dx*0.5d0,dy*0.5d0,dz*0.5d0
      write(10,"('ASPECT_RATIO ',e11.4,' ',e11.4,' ',e11.4)")dx,dy,dz
      write(10,"('')")
      write(10,"('POINT_DATA ',I0)")mi*mj*mk
      write(10,"('SCALARS p float')")
      write(10,"('LOOKUP_TABLE default')")
      do k=1,mk
      do j=1,mj
      do i=1,mi
      vorx00=(
     1  vorx(i,j-1,k-1)
     2 +vorx(i,j  ,k-1)
     3 +vorx(i,j-1,k  )
     4 +vorx(i,j  ,k  ))*0.25d0
      q00=q(i,j,k)
      if(vorx00.ge.0.0d0)q00=0.0d0
      write(10,"(e10.3)")q00
      enddo
      enddo
      enddo
      close(10)
      write(*,*)'output ',fname

      return
      end
