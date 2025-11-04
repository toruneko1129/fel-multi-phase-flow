      subroutine mkvtk_velocity_nj(svall, nstep, dx, dy, dz, u, v, w)
      implicit none
      integer svall(3), nstep
      real*8  dx, dy, dz
      ! スタガード速度（ゴースト付でOK）
      real*8  u(-2:svall(1)+3, -2:svall(2)+3, -2:svall(3)+3)
      real*8  v(-2:svall(1)+3, -2:svall(2)+3, -2:svall(3)+3)
      real*8  w(-2:svall(1)+3, -2:svall(2)+3, -2:svall(3)+3)

      integer :: ni, nj, nk, i, j, k
      real*8  :: uc, vc, wc
      character*32 :: fname
      integer :: ncell

      ni = svall(1)   ! セル数（x方向）
      nj = svall(2)   ! セル数（y方向）
      nk = svall(3)   ! セル数（z方向）
      ncell = ni*1*nk

      write(fname,'("vel(nj)_",i7.7,".vtk")') nstep
      open(10, file=fname, status='replace')

      ! ---- VTK Legacy / STRUCTURED_POINTS / ASCII ----
      write(10,'("# vtk DataFile Version 2.0")')
      write(10,'("velocity at cell centers, step=",I0)') nstep
      write(10,'("ASCII")')
      write(10,'("DATASET STRUCTURED_POINTS")')
      write(10,'("DIMENSIONS ",I0," ",I0," ",I0)') ni+1, 2, nk+1
      write(10,'("ORIGIN ",ES22.14," ",ES22.14," ",ES22.14)') 0.0d0, 0.0d0, 0.0d0
      ! 現行仕様では SPACING を使うのが推奨
      write(10,'("SPACING ",ES22.14," ",ES22.14," ",ES22.14)') dx, dy, dz
      write(10,'("")')

      ! ---- セルデータ（コロケート済みの速度ベクトル）----
      write(10,'("CELL_DATA ",I0)') ncell
      write(10,'("VECTORS velocity float")')

      do k = 1, nk
        do j = nj, nj
          do i = 1, ni
            ! スタガード → セル中心（コロケート）
            uc = 0.5d0 * ( u(i-1,j  ,k  ) + u(i  ,j  ,k  ) )
            vc = 0.5d0 * ( v(i  ,j-1,k  ) + v(i  ,j  ,k  ) )
            wc = 0.5d0 * ( w(i  ,j  ,k-1) + w(i  ,j  ,k  ) )
            write(10,'(3ES22.14)') uc, vc, wc
          end do
        end do
      end do

      close(10)
      write(*,*) 'output ', trim(fname)
      return
      end
