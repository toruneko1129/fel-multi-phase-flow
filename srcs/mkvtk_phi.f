      subroutine mkvtk_phi(svall,nstep,dx,dy,dz,q)
      implicit none
      integer svall(3),nstep
      real*8 dx,dy,dz
      real*8 q(-2:svall(1)+3,-2:svall(2)+3,-2:svall(3)+3)

      character*32 fname
      integer i,j,k,mi,mj,mk
      real*8 q00

      mi=svall(1)
      mj=svall(2)
      mk=svall(3)

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(j,k)
!$OMP$ SHARED(mi,mj,mk,q)      
      do k=0,mk+1
      do j=0,mj+1
      q(   0,j,k)=q(mi,j,k)
      q(mi+1,j,k)=q( 1,j,k)
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,k)
!$OMP$ SHARED(mi,mj,mk,q)      
      do k=0,mk+1
      do i=0,mi+1
      q(i,   0,k)=0.0d0
      q(i,mj+1,k)=0.0d0
!      q(i,   0,k)=q(i, 1,k)
!      q(i,mj+1,k)=q(i,mj,k)
      enddo
      enddo
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j)
!$OMP$ SHARED(mi,mj,mk,q)      
      do j=0,mj+1
      do i=0,mi+1
!      q(i,j,   0)=1.0d0
!      q(i,j,mk+1)=1.0d0
      q(i,j,   0)=q(i,j, 1)
      q(i,j,mk+1)=q(i,j,mk)
      enddo
      enddo
!$OMP  END PARALLEL DO

      write(fname,'("phi_",i7.7,".vtk")')nstep
      open(10,file=fname,status='replace')
      write(10,"('# vtk DataFile Version 2.0')")
      write(10,"('t=',I0)")nstep
      write(10,"('ASCII')")
      write(10,"('DATASET STRUCTURED_POINTS')")
      write(10,"('DIMENSIONS ',I0,' ',I0,' ',I0)")mi+1,mj+1,mk+1
      write(10,"('ORIGIN ',e11.4,' ',e11.4,' ',e11.4)")0.0,0.0,0.0
      write(10,"('ASPECT_RATIO ',e11.4,' ',e11.4,' ',e11.4)")dx,dy,dz
      write(10,"('')")
      write(10,"('POINT_DATA ',I0)")(mi+1)*(mj+1)*(mk+1)
      write(10,"('SCALARS p float')")
      write(10,"('LOOKUP_TABLE default')")
      do k=0,mk
      do j=0,mj
      do i=0,mi
      q00=(
     1 +q(i  ,j  ,k  )
     2 +q(i+1,j  ,k  )
     3 +q(i  ,j+1,k  )
     4 +q(i+1,j+1,k  )
     5 +q(i  ,j  ,k+1)
     6 +q(i+1,j  ,k+1)
     7 +q(i  ,j+1,k+1)
     8 +q(i+1,j+1,k+1) )*0.125d0
      write(10,"(e10.3)")q00
      enddo
      enddo
      enddo
      close(10)
      write(*,*)'output ',fname

      return
      end
