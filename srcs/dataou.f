      subroutine dataou(ipara,ID,ni,nj,nk,nbub,nstep,time
     & ,u,v,w,p,uo,vo,wo,po,phi)

      implicit none
      integer ipara,ID
      integer ni,nj,nk,nbub,nstep
      real*8 time
      real*8    u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    p(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   uo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   vo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   wo(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   po(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  phi(-2:ni+3,-2:nj+3,-2:nk+3,0:nbub)
      integer i,j,k,l
      character*32 fname

      if(ipara.eq.1)then
      write(fname,'("x",i7.7,"_",i4.4)')nstep,ID
      else
      write(fname,'("x",i7.7)')nstep
      endif

      open(10,file=fname,form='unformatted')
      write(10)
     1   (((   u(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     2 , (((   v(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     3 , (((   w(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     4 , (((   p(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     5 , (((  uo(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     6 , (((  vo(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     7 , (((  wo(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     8 , (((  po(i,j,k),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3)
     9 , (((( phi(i,j,k,l),i=-2,ni+3),j=-2,nj+3),k=-2,nk+3),l=0,nbub)
     & ,time
     1 ,nstep
      close(10)
      write(*,*)'output ',fname
c
      return
      end
