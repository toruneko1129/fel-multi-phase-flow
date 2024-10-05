      write(fname,'("uwphi",i7.7)')nstep
      open(10,file=fname)
      do k=0,svall(3)
      do i=0,svall(1)
      x00=dble(i)*dx
      z00=dble(k)*dz
      u00=  u_slice(i,k)
      w00=  w_slice(i,k)
      q00=phi_slice(i,k)
      if(abs(u00).lt.1.0d-40)u00=0.0d0
      if(abs(w00).lt.1.0d-40)w00=0.0d0
      if(abs(q00).lt.1.0d-40)q00=0.0d0
      write(10,'(20e20.10)') &
       x00 &
      ,z00 &
      ,u00 &
      ,w00 &
      ,q00
      enddo
      write(10,'()')
      enddo
      close(10)
      write(*,*)'output ',fname
