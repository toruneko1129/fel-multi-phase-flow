        subroutine mk_vel(nstep,time,vel,re,velocity,l)
        implicit none
        integer nstep,l
        integer, parameter:: n = 3
        real*8 dt,vel(n),re(n),velocity

        integer, parameter:: funit = 71
        integer i
        real*8 time
        character*32  filename
        write(filename,'("velocity_",i2.2".dat")')l
       
       open(unit=funit,file=filename
     & ,status="unknown",position="append")
       write(funit,'(20e20.10)')time
     & ,vel(1),vel(2),vel(3),re(1),re(2),velocity
       close(funit)
       
       return
       end
