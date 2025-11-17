      subroutine mkvtk_phi_cell(svall, nstep, dx, dy, dz, q)
!----- セル中心の φ を VTK Legacy (STRUCTURED_POINTS) の CELL_DATA で出力
!      svall(1:3) = (mi, mj, mk)  … セル数
!      配列 q はゴースト付きだが、出力は内点 i=1..mi, j=1..mj, k=1..mk のみ
      implicit none
      integer           svall(3), nstep
      real*8            dx, dy, dz
      real*8            q(-2:svall(1)+3, -2:svall(2)+3, -2:svall(3)+3)

      integer           mi, mj, mk
      integer           i, j, k
      character*64      fname

      mi = svall(1)
      mj = svall(2)
      mk = svall(3)

      write(fname,'("phi_",i7.7,".vtk")') nstep
      open(10, file=fname, status='replace', action='write')

!---- VTK Legacy ヘッダ
      write(10,'(A)') '# vtk DataFile Version 2.0'
      write(10,'(A)') 'phi (cell data)'
      write(10,'(A)') 'ASCII'
      write(10,'(A)') 'DATASET STRUCTURED_POINTS'
!     DIMENSIONS は「点数」= (mi+1, mj+1, mk+1)
      write(10,'("DIMENSIONS ",I0," ",I0," ",I0)') mi+1, mj+1, mk+1
!     原点と格子ピッチを明示（ASPECT_RATIO ではなく SPACING を推奨）
      write(10,'("ORIGIN ",ES16.9," ",ES16.9," ",ES16.9)') 0.0d0, 0.0d0, 0.0d0
      write(10,'("SPACING ",ES16.9," ",ES16.9," ",ES16.9)') dx, dy, dz
      write(10,'(A)') ''

!---- CELL_DATA ブロック（件数は mi*mj*mk）
      write(10,'("CELL_DATA ",I0)') mi*mj*mk
      write(10,'(A)') 'SCALARS phi float 1'
      write(10,'(A)') 'LOOKUP_TABLE default'

!---- データ本体（VTKは x→y→z の順で i が最内ループ）
      do k = 1, mk
        do j = 1, mj
          do i = 1, mi
            write(10,'(ES16.9)') q(i,j,k)
          enddo
        enddo
      enddo

      close(10)
      write(*,'(A)') 'output '//trim(fname)
      return
      end
