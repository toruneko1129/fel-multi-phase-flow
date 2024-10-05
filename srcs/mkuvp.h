      write(fname,'("uvwp",i7.7)')nstep
      open(10,file=fname)
      do j=0,svall(2)
      do i=0,svall(1)
      x00=dble(i)*dx
      y00=dble(j)*dy
      u00=  u_slice(i,j)
      v00=  v_slice(i,j)
      q00=phi_slice(i,j)
      if(abs(u00).lt.1.0d-40)u00=0.0d0
      if(abs(v00).lt.1.0d-40)v00=0.0d0
      if(abs(q00).lt.1.0d-40)q00=0.0d0
      write(10,'(20e20.10)') &
       x00 &
      ,y00 &
      ,0.0d0 &
      ,u00 &
      ,v00 &
      ,q00
      enddo
      write(10,'()')
      enddo
      close(10)
      write(*,*)'output ',fname
