        subroutine mk_center(nstep,time,center,velocity,l)
        implicit none
        integer nstep,l
        integer, parameter:: n = 3
        real*8 dt,center(n)
        real*8 vector(n),velocity

        integer, parameter:: funit = 71
        integer i
        real*8 time
        character*32 filename

        write(filename,'("center_",i2.2,".dat")')l

       open(unit=funit,file=filename
     & ,status="unknown",position="append")
       write(funit,'(20e20.10)')time
     & ,center(1),center(2),center(3),velocity
     & ,sqrt((center(1)-0.5d0)**2+(center(3)-0.5d0)**2) 
       close(funit)
       
       return
       end
