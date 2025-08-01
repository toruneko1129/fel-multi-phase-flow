program main
implicit none
include 'mpif.h'
include 'param.h'

!nsv: number of mesh
!svall: use svn to determine mesh resolution
integer :: nsv
integer :: svall(3)
integer :: ierr,ndiv,ID,ipara
integer :: nID(6)
integer,dimension(:),allocatable :: key
real(8),dimension(:),allocatable :: sendbuf,recvbuf
real(8),dimension(:),allocatable :: sendjb,recvjb
real(8),dimension(:),allocatable :: dyinv_array

integer irestart
integer ni,nj,nk
integer nbub
parameter (nbub=1)              !change 
real(8) :: xxc(nbub),yyc(nbub),zzc(nbub)

real(8) :: pi,cfl,time,dt,xl,yl,zl,dx,dy,dz,dxinv,dyinv,dzinv
real(8) :: surface_tension,rhol,rhog,rmul,rmug,grv,grvb,grvp,angle_deg,angle_rad,uwall
real(8) :: l1_a,l2_a,theta_0_a,zeta_a
real(8) :: l1_b,l2_b,theta_0_b,zeta_b
real(8) :: bet_mthinc
real(8) :: particle_radius,particle_init_x,particle_init_y,particle_init_z
integer nmax,idout,imkuvp,imkvtk,ibudget,imon_t,nstep,nstep0,tscale

real(8),dimension(:,:,:),allocatable :: u,v,w,uo,vo,wo,un,vn,wn,p,po,pn,phat,dp,div
real(8),dimension(:,:,:),allocatable :: rho,rmu,rhon,rmun
real(8),dimension(:,:,:),allocatable :: adv_u ,adv_v ,adv_w
real(8),dimension(:,:,:),allocatable :: adv_uo,adv_vo,adv_wo
real(8),dimension(:,:,:),allocatable :: prs_u ,prs_v ,prs_w
real(8),dimension(:,:,:),allocatable :: vis_u ,vis_v ,vis_w
real(8),dimension(:,:,:),allocatable :: vis_un,vis_vn,vis_wn
real(8),dimension(:,:,:),allocatable :: fst_u ,fst_v ,fst_w
real(8),dimension(:,:,:),allocatable :: fst_un,fst_vn,fst_wn
real(8),dimension(:,:,:),allocatable :: sum_fst_u ,sum_fst_v ,sum_fst_w
real(8),dimension(:,:,:),allocatable :: sum_fst_un,sum_fst_vn,sum_fst_wn
real(8),dimension(:,:,:),allocatable :: src_u,src_v,src_w
real(8),dimension(:,:,:),allocatable :: phix,phiy,phiz,flphix,flphiy,flphiz
real(8),dimension(:,:,:,:),allocatable :: phi,phin,phil_all
real(8),dimension(:,:,:,:),allocatable :: s ,tau,sn,taun
real(8),dimension(:,:,:,:),allocatable ::  rn_pp, rn_12, rn_13, rn_23,rkap
real(8),dimension(:,:,:,:),allocatable :: rnn_pp,rnn_12,rnn_13,rnn_23,rkapn

real(8),dimension(:,:,:),allocatable :: au_w_u,au_e_u,au_s_u,au_n_u,au_b_u,au_t_u,au_p_u
real(8),dimension(:,:,:),allocatable :: av_ws_u,av_es_u,av_wn_u,av_en_u,aw_wb_u,aw_eb_u,aw_wt_u,aw_et_u
real(8),dimension(:,:,:),allocatable :: av_w_v,av_e_v,av_s_v,av_n_v,av_b_v,av_t_v,av_p_v
real(8),dimension(:,:,:),allocatable :: au_sw_v,au_nw_v,au_se_v,au_ne_v,aw_sb_v,aw_nb_v,aw_st_v,aw_nt_v
real(8),dimension(:,:,:),allocatable :: aw_w_w,aw_e_w,aw_s_w,aw_n_w,aw_b_w,aw_t_w,aw_p_w
real(8),dimension(:,:,:),allocatable :: au_bw_w,au_tw_w,au_be_w,au_te_w,av_bs_w,av_ts_w,av_bn_w,av_tn_w
real(8),dimension(:,:,:),allocatable :: aw_p,ae_p,as_p,an_p,ab_p,at_p,ap_p

real(8),dimension(:),allocatable :: fftdata,atdma,btdma_r,btdma_i
real(8),dimension(:),allocatable :: as_p_fft,an_p_fft,ap_p_fft
real(8),dimension(:),allocatable :: as_u_fft,an_u_fft,ap_u_fft
real(8),dimension(:),allocatable :: as_v_fft,an_v_fft,ap_v_fft
real(8),dimension(:),allocatable :: ab_w_fft,at_w_fft,ap_w_fft
real(8),dimension(:,:,:),allocatable :: src_wkdiv,wnkdiv

real(8),dimension(:,:,:),allocatable :: phir_r,phir_i
real(8),dimension(:,:,:),allocatable :: phiw_r,phiw_i

real(8),dimension(:,:,:),allocatable :: theta_0_array
real(8),dimension(:,:,:),allocatable :: l1_array
real(8),dimension(:,:,:),allocatable :: l2_array
real(8),dimension(:,:,:),allocatable :: zeta_array
real(8),dimension(:,:,:),allocatable :: theta_array
integer period, ratio_a

real(8),dimension(:,:,:),allocatable :: vorx,q
real(8),dimension(:,:,:),allocatable :: phi_all,q_all,vorx_all,wrk_all,p_all
real(8),dimension(:,:),allocatable :: u_slice,v_slice,w_slice,phi_slice

real(8),dimension(:),allocatable :: center,vel,re

integer i,j,k,l,ii,jj,kk,area_pre1,area_pre2,lap
character(32) fname
real(8) :: x00,y00,z00,u00,v00,w00,p00,q00
real(8) :: avphin,err_div0,err_div
real(8) :: rho_av,rhon_av
real(8) :: center_pre1,center_pre2,velocity

! the number of grid points over the entire region
! in Legendre case, use (8svn, svn, 2)
!!!!
nsv=128

svall(1)=nsv*8
svall(2)=nsv
svall(3)=2

! irestart=1 if computation will be restarted (input data are needed).
! irestart=0 if computation will be started from t=0.

irestart=0

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,ndiv,ierr)
call mpi_comm_rank(mpi_comm_world,ID,ierr)
ipara=1
if(ndiv.eq.1)ipara=0

if(mod(svall(1),ndiv).ne.0.or.mod(svall(2),ndiv).ne.0)then
  if(ID.eq.0)then
    print *,'Invalid number of PE'
    print *,'set PE (=ndiv) to satisfy svall(1)%ndiv=0 and svall(2)%ndiv=0'
  endif
  call mpi_finalize(ierr)
  stop
endif

if(svall(1).ne.2**int(log(dble(svall(1)+1))/log(2.0d0)).or. &
   svall(2).ne.2**int(log(dble(svall(2)+1))/log(2.0d0)).or. &
   svall(3).ne.2**int(log(dble(svall(3)+1))/log(2.0d0)))then
  if(ID.eq.0)then
    print *,'Invalid number of grid points'
    print *,'svall(m) for each m (=1, 2, or 3) must be 2**n with integer n'
  endif
  call mpi_finalize(ierr)
  stop
endif

nID(X_MINUS)=-1
nID(X_PLUS )=-1
nID(Y_MINUS)=-1
nID(Y_PLUS )=-1
nID(Z_MINUS)=-1
nID(Z_PLUS )=-1
if(ID.ge.1     )nID(Y_MINUS)=ID-1
if(ID.ne.ndiv-1)nID(Y_PLUS )=ID+1

ni=svall(1)
nj=svall(2)/ndiv
nk=svall(3)

!open(10,file='inp')
!read(10,*)scaled_init_y
!close(10)

pi=atan(1.0d0)*4.0d0

! drv_prs_grad: driving pressure gradient in x
! xl, yl, zl: lengthes in x, y and z directions to describe domain size
! surface_tension: surface_tension
! particle_radius: particle radius
! particle_init_[xyz]: initial position from the origin at the system centroid
! ni, nj, nk: the numbers of grid points
! dx, dy, dz: grid widths
! uwall: wall velocity
! theta_0: static contact angle at the wall[deg]

!!!!
xl=1.36d2
yl=1.36d1
zl=4.25

rhol=8.1d-1
rhog=8.1d-1
rmul=1.95d0
rmug=1.95d0
surface_tension=5.5d0

uwall = 0.25d0
l1_a = 2.165d0
l2_a = l1_a
theta_0_a = 90.0d0
zeta_a = 0.21d0 * 6.0d0 * (rmul / l1_a + rmug / l2_a)

l1_b = 2.379d0
l2_b = l1_b * 3.67d0/1.625d0
theta_0_b = 64.0d0
zeta_b = 0.21d0 * 6.0d0 * (rmul / l1_b + rmug / l2_b)

!pattern width
period = 8
ratio_a = 0

!calculation gravity 
!!!!
grv =0.0d0
angle_deg=0.0d0
angle_rad=angle_deg*pi/180.0d0
grvb=grv*sin(angle_rad)
grvp=-grv*cos(angle_rad)
if(ID.eq.0)then
  write(*,'("deg,rad,grvb,grvp",20e20.10)')angle_deg,angle_rad,grvb,grvp
endif

dx=xl/dble(svall(1))
dy=yl/dble(svall(2))
dz=zl/dble(svall(3))

! allocate memory for each array variable
include'allocate.h'

!init static contact angle at the wall
call init_array_pt_ratio(ni,nj,nk,theta_0_a,theta_0_b,theta_0_array,period,ratio_a)
call init_array_pt_ratio(ni,nj,nk,l1_a,l1_b,l1_array,period,ratio_a)
call init_array_pt_ratio(ni,nj,nk,l2_a,l2_b,l2_array,period,ratio_a)
call init_array_pt_ratio(ni,nj,nk,zeta_a,zeta_b,zeta_array,period,ratio_a)

dxinv=1.0d0/dx
dyinv=1.0d0/dy

dyinv_array(-2)=1.0d0/dy
dyinv_array(-1)=1.0d0/dy
dyinv_array(0)=1.0d0/dy
dyinv_array(1)=1.0d0/dy
dyinv_array(svall(2)+1)=1.0d0/dy
dyinv_array(svall(2)+2)=1.0d0/dy
dyinv_array(svall(2)+3)=1.0d0/dy
do i=2,svall(2)
dyinv_array(i)=1.0d0/dy
enddo

dzinv=1.0d0/dz


!cfl: CFL number

cfl=0.05d0
bet_mthinc=2.0d0
!bet_ibm=2.0d0

!ccc
!ccc nmax    the maximum number of computational time steps
!ccc idout   interval for outputing the backup data
!ccc imkuvp  interval for outputing the instantaneous field
!ccc imon_t  interval for monitoring the time step
!ccc ibudget interval for writing budgets
!ccc

tscale  =1.0d0
nmax    =12000*nsv/32/tscale
idout   =1200000
imkuvp  =1000000
imkvtk  =nmax/120
imon_t  =nmax/120
ibudget =imon_t

time=0.0d0

if(ID.eq.0)then
write(*,'("svall(1)=           ",1i9)')svall(1)
write(*,'("svall(2)=           ",1i9)')svall(2)
write(*,'("svall(3)=           ",1i9)')svall(3)
write(*,'("ndiv=               ",1i9)')ndiv
write(*,'("ni=                 ",1i9)')ni
write(*,'("nj=                 ",1i9)')nj
write(*,'("nk=                 ",1i9)')nk
write(*,'("irestart=           ",1i9)')irestart
write(*,*)
write(*,'("xl=                 ",20e20.10)')xl
write(*,'("yl=                 ",20e20.10)')yl
write(*,'("zl=                 ",20e20.10)')zl
write(*,'("cfl                 ",20e20.10)')cfl
write(*,'("surface_tension     ",20e20.10)')surface_tension
write(*,'("rhol                ",20e20.10)')rhol
write(*,'("rhog                ",20e20.10)')rhog
write(*,'("rmul                ",20e20.10)')rmul
write(*,'("rmug                ",20e20.10)')rmug
write(*,'("grv                 ",20e20.10)')grv
write(*,'("bet_mthinc          ",20e20.10)')bet_mthinc
write(*,'("uwall               ",20e20.10)')uwall
write(*,'("l1_a                ",20e20.10)')l1_a
write(*,'("l2_a                ",20e20.10)')l2_a
write(*,'("theta_0_a           ",20e20.10)')theta_0_a
write(*,'("zeta_a              ",20e20.10)')zeta_a
write(*,'("l1_b                ",20e20.10)')l1_b
write(*,'("l2_b                ",20e20.10)')l2_b
write(*,'("theta_0_b           ",20e20.10)')theta_0_b
write(*,'("zeta_b              ",20e20.10)')zeta_b
write(*,'("pattern period      ",1i9)')period
write(*,'("pattern ratio_a     ",1i9)')ratio_a
write(*,'("pattern ratio_b     ",1i9)')period-ratio_a
write(*,*)
write(*,'("nmax                ",1i9)')nmax
write(*,'("idout               ",1i9)')idout
write(*,'("imkuvp              ",1i9)')imkuvp
write(*,'("ibudget             ",1i9)')ibudget
write(*,'("imon_t              ",1i9)')imon_t
write(*,*)
endif

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)

call init(ni,nj,nk,u,v,w,p,uo,vo,wo,po,dp,phi,nbub)
call cpy(ni,nj,nk,u,un)
call cpy(ni,nj,nk,v,vn)
call cpy(ni,nj,nk,w,wn)

particle_radius=0.5d-3
!particle_radius=2.5d-3

if(ID.eq.0)then
write(*,'("xl,yl,zl",10e20.10)')xl,yl,zl
endif

if(irestart.eq.0)then
write(*,'("START")')

!for tow phase flow, nbub must be 1
call init_phi(ndiv, svall, phi, nbub, 1, nsv)


do l=1,nbub
!>contact angle condition
!call bnd_neumann(nID,ni,nj,nk,phi(-2,-2,-2,l))
call gnbc(nID, ni, nj, nk, u, uwall, theta_0_array, &
          surface_tension, zeta_array, theta_array)
call bnd_contact_angle(nID,ni,nj,nk,phi(-2,-2,-2,l),theta_array,dx,dy)
!call bnd_dirichlet(nID,ni,nj,nk,phi(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,phi(-2,-2,-2,l))
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,phi(-2,-2,-2,l))
enddo

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)
call summation(ni,nj,nk,phi,nbub)

call bndu(nID,ni,nj,nk,u ,v ,w ,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
call bndu(nID,ni,nj,nk,uo,vo,wo,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,u )
call bnd_periodic(ni,nj,nk,v )
call bnd_periodic(ni,nj,nk,w )
call bnd_periodic(ni,nj,nk,uo)
call bnd_periodic(ni,nj,nk,vo)
call bnd_periodic(ni,nj,nk,wo)
call bnd_periodic(ni,nj,nk,un)
call bnd_periodic(ni,nj,nk,vn)
call bnd_periodic(ni,nj,nk,wn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,u )
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,v )
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,w )
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,uo)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vo)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wo)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,un)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wn)

!>contact angle condition
!call bnd_neumann(nID,ni,nj,nk,phi(-2,-2,-2,0))
call gnbc(nID, ni, nj, nk, u, uwall, theta_0_array, &
          surface_tension, zeta_array, theta_array)
call bnd_contact_angle(nID,ni,nj,nk,phi(-2,-2,-2,0),theta_array,dx,dy)
!call bnd_dirichlet(nID,ni,nj,nk,phi(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,phi(-2,-2,-2,0))
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,phi(-2,-2,-2,0))

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)
endif
!nomal start end

nstep=0

!cc<  restart
if(irestart.eq.1)then
  write(*,'("RESTART")')
  call datain(ipara,ID,ni,nj,nk,nbub,nstep,time,u,v,w,p,uo,vo,wo,po,phi)
  call bndu(nID,ni,nj,nk,u ,v ,w ,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
  call bndu(nID,ni,nj,nk,uo,vo,wo,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
  call bnd_periodic(ni,nj,nk,u )
  call bnd_periodic(ni,nj,nk,v )
  call bnd_periodic(ni,nj,nk,w )
  call bnd_periodic(ni,nj,nk,p )
  call bnd_periodic(ni,nj,nk,uo)
  call bnd_periodic(ni,nj,nk,vo)
  call bnd_periodic(ni,nj,nk,wo)
  call bnd_periodic(ni,nj,nk,po)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,u )
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,v )
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,w )
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,p )
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,uo)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vo)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wo)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,po)
  !>contact angle condition
  !call bnd_neumann(nID,ni,nj,nk,phi )
  call gnbc(nID, ni, nj, nk, u, uwall, theta_0_array, &
            surface_tension, zeta_array, theta_array)
  call bnd_contact_angle(nID,ni,nj,nk,phi,theta_array,dx,dy)
  !call bnd_dirichlet(nID,ni,nj,nk,phi)
  call bnd_periodic(ni,nj,nk,phi )
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,phi )
endif
!cc< restart end

j = particle_radius/dy
call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,0,1,1,wrk_all,  u,  u_slice)
call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,1,1,0,wrk_all,  w,  w_slice)
call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,1,1,1,wrk_all,phi,phi_slice)
if(ID.eq.0)then
  include'mkuwp.h'
endif

call calvorx(ni,nj,nk,dyinv,dzinv,v,w,vorx)
call bnd_periodic(ni,nj,nk,vorx)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vorx)
call cal_advu(ni,nj,nk,dxinv,dyinv,dzinv,u,v,w,adv_u,adv_v,adv_w)
call calq(ni,nj,nk,dxinv,dyinv,dzinv,adv_u,adv_v,adv_w,q)
call bnd_periodic(ni,nj,nk,q)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,q)

do l=1,nbub
call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,phi(-2,-2,-2,l),phil_all(-2,-2,-2,l))
enddo
call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all, phi, phi_all)
call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,   q,   q_all)
call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,vorx,vorx_all)
call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,   p,   p_all)
if(irestart.eq.0)then
  vel(1)= 0.0d0
  vel(2)= 0.0d0
  vel(3)= 0.0d0
  re(1) = 0.0d0
  re(2) = 0.0d0
  velocity = 0.0d0
endif
if(ID.eq.0)then
  !do l=1,nbub
  !call mkvtk_phil(svall,nstep,dx,dy,dz,phil_all(-2,-2,-2,l),l)
  !enddo
  call mkvtk_phi(svall,nstep,dx,dy,dz,phi_all)
  call mkvtk_p(svall,nstep,dx,dy,dz,p_all)
  !call   mkvtk_q(svall,nstep,dx,dy,dz,vorx_all,q_all)
  do l=1,nbub
  center(1)= xxc(l)
  center(2)= yyc(l)
  center(3)= zzc(l)
  call mk_vel(nstep,time,vel,re,velocity,l)
  call mk_center(nstep,time,center,velocity,l)
  enddo
endif
call mpi_barrier(mpi_comm_world,ierr)
call flush(6)
!stop

nstep0=nstep

if(ID.eq.0)then
  write(*,'()')
  write(*,'("nstep0     ",1i9)')nstep0
  write(*,'("time0      ",20e20.10)')time
  write(*,'()')
endif

call cal_coef_prs(svall(2),rhog,dxinv,dyinv,dzinv,aw_p,ae_p,as_p,an_p,ab_p,at_p,ap_p)
call cal_coef_solp_fft(svall(2),as_p,an_p,as_p_fft,an_p_fft,ap_p_fft)

!ccc
!ccc<main routine
!ccc
do nstep=nstep0+1,nstep0+nmax
if(mod(nstep,imon_t).eq.0.and.ID.eq.0)then
write(*,*)'---------------------------------------'
write(*,'("nstep= ",1i9.9)')nstep
endif

!call caldt(ipara,nID,ID,ndiv,ni,nj,nk,nstep,imon_t,dxinv,dyinv,dzinv,cfl,rhol,rhog,rmul,rmug,surface_tension,u,v,w,dt,time)
!>tmp changed
dt=32.0d-2/nsv*tscale
time=time+dt
call mpi_barrier(mpi_comm_world,ierr)
if(mod(nstep,imon_t).eq.0.and.ID.eq.0)then
write(*,'("time=",1e17.10," dt=",1e17.10)'), time, dt
endif
call flush(6)

do l=1,nbub
!cc< mthinc
call cal_grad_p2a(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phi(-2,-2,-2,l),phix,phiy,phiz)
call solphi_mthinc1(ni,nj,nk,dxinv,dyinv,dzinv,dt,bet_mthinc,phi(-2,-2,-2,l),phix,phiy,phiz,u,v,w,flphix,flphiy,flphiz)
call solphi_mthinc2(ni,nj,nk,dxinv,dyinv,dzinv,dt,u,v,w,flphix,flphiy,flphiz,phi(-2,-2,-2,l),phin(-2,-2,-2,l))
call cal_grad_p2a(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phin(-2,-2,-2,l),phix,phiy,phiz)
call solphi_mthinc3(ipara,ni,nj,nk,dxinv,dyinv,dzinv,bet_mthinc,phix,phiy,phiz,phi(-2,-2,-2,l),phin(-2,-2,-2,l))
!>contact angle condition
!call bnd_neumann(nID,ni,nj,nk,phin(-2,-2,-2,l))
call gnbc(nID, ni, nj, nk, u, uwall, theta_0_array, &
          surface_tension, zeta_array, theta_array)
call bnd_contact_angle(nID,ni,nj,nk,phin(-2,-2,-2,l),theta_array,dx,dy)
!call bnd_dirichlet(nID,ni,nj,nk,phin(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,phin(-2,-2,-2,l))
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,phin(-2,-2,-2,l))
enddo

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)

call summation(ni,nj,nk,phin,nbub)
!>contact angle condition
call bnd_neumann(nID,ni,nj,nk,phin(-2,-2,-2,0))
call gnbc(nID, ni, nj, nk, u, uwall, theta_0_array, &
          surface_tension, zeta_array, theta_array)
call bnd_contact_angle(nID,ni,nj,nk,phi(-2,-2,-2,0),theta_array,dx,dy)
!call bnd_dirichlet(nID,ni,nj,nk,phi(-2,-2,-2,0))
call bnd_periodic(ni,nj,nk,phin(-2,-2,-2,0))
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,phin(-2,-2,-2,0))

if(mod(nstep,imon_t).eq.0)then
  call cal_av(ipara,ndiv,ni,nj,nk,phin,avphin)
  if(ID.eq.0)then
    write(*,'("AV_PHI   ",1i10,20e20.10)')nstep,avphin
  endif
endif

call cal_arith_mean(ni,nj,nk,phi ,rhol,rhog,rho )
call cal_arith_mean(ni,nj,nk,phin,rhol,rhog,rhon)
call cal_arith_mean(ni,nj,nk,phi ,rmul,rmug,rmu )
call cal_arith_mean(ni,nj,nk,phin,rmul,rmug,rmun)

!cc< WENO
call cal_advu_weno5(ni,nj,nk,dxinv,dyinv,dzinv,u ,v ,w ,adv_u ,adv_v ,adv_w )
call cal_advu_weno5(ni,nj,nk,dxinv,dyinv,dzinv,uo,vo,wo,adv_uo,adv_vo,adv_wo)

call calphat(ni,nj,nk,p,po,phat)
call cal_grad_prs(ni,nj,nk,dxinv,dyinv,dzinv,phat,prs_u,prs_v,prs_w)

call calsij(ni,nj,nk,dxinv,dyinv,dzinv,u ,v ,w ,s )
call cal_arith_tau(ni,nj,nk,s,rmu,tau)

call cal_div_tensor(ni,nj,nk,dxinv,dyinv,dzinv,tau ,vis_u ,vis_v ,vis_w )


!calculate surface tension force (sum_fst_[uvw],sumf_fst_[uvw]n)
call init_q(ni,nj,nk,sum_fst_u,sum_fst_v,sum_fst_w)
call init_q(ni,nj,nk,sum_fst_un,sum_fst_vn,sum_fst_wn)
do l=1,nbub
call cal_grad_p2a(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phi(-2,-2,-2,l),phix,phiy,phiz)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,1,1,1,phix,phiy,phiz,rn_pp)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,0,0,1,phix,phiy,phiz,rn_12)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,0,1,0,phix,phiy,phiz,rn_13)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,1,0,0,phix,phiy,phiz,rn_23)

call cal_grad_p2a(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phin(-2,-2,-2,l),phix,phiy,phiz)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,1,1,1,phix,phiy,phiz,rnn_pp)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,0,0,1,phix,phiy,phiz,rnn_12)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,0,1,0,phix,phiy,phiz,rnn_13)
call calrn(ni,nj,nk,dxinv,dyinv,dzinv,1,0,0,phix,phiy,phiz,rnn_23)

call calrkap(ni,nj,nk,dxinv,dyinv,dzinv, rn_pp, rn_12, rn_13, rn_23,rkap )
call calrkap(ni,nj,nk,dxinv,dyinv,dzinv,rnn_pp,rnn_12,rnn_13,rnn_23,rkapn)

call cal_grad_p2uvw(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phi(-2,-2,-2,l),phix,phiy,phiz)
call calfst(ni,nj,nk,rhol,rhog,surface_tension,phix,phiy,phiz,rkap,fst_u,fst_v,fst_w)

call cal_grad_p2uvw(ID,svall(2),ni,nj,nk,dxinv,dyinv_array,dzinv,phin(-2,-2,-2,l),phix,phiy,phiz)
call calfst(ni,nj,nk,rhol,rhog,surface_tension,phix,phiy,phiz,rkapn,fst_un,fst_vn,fst_wn)
call sum_fst(ni,nj,nk,fst_u,fst_v,fst_w,sum_fst_u,sum_fst_v,sum_fst_w)
call sum_fst(ni,nj,nk,fst_un,fst_vn,fst_wn,sum_fst_un,sum_fst_vn,sum_fst_wn)
enddo
call cal_av(ipara,ndiv,ni,nj,nk,rho ,rho_av)
call cal_av(ipara,ndiv,ni,nj,nk,rhon,rhon_av)
rho_av=(rho_av+rhon_av)*0.5d0

!if(ID.eq.0)then
!     write(*,'("sum_fst",10e20.10)')sum_fst_u(24,0,32),sum_fst_u(24,1,32),sum_fst_u(24,nj,32)
!endif

call cal_srcu(nID,ni,nj,nk &
  ,rho_av,dt,grvb,grvp &
  ,u,v,w,rho,rhon  &
  ,adv_u ,adv_v ,adv_w  &
  ,adv_uo,adv_vo,adv_wo &
  ,prs_u ,prs_v ,prs_w  &
  ,vis_u ,vis_v ,vis_w  &
  ,sum_fst_u ,sum_fst_v ,sum_fst_w  &
  ,sum_fst_un,sum_fst_vn,sum_fst_wn &
  ,src_u ,src_v ,src_w )

call cal_arith_coef_vis(ni,nj,nk,dxinv,dyinv,dzinv &
 ,    rmun &
 , au_w_u, au_e_u, au_s_u, au_n_u &
 , au_b_u, au_t_u, au_p_u         &
 ,av_ws_u,av_es_u,av_wn_u,av_en_u &
 ,aw_wb_u,aw_eb_u,aw_wt_u,aw_et_u &
 , av_w_v, av_e_v, av_s_v, av_n_v &
 , av_b_v, av_t_v, av_p_v         &
 ,au_sw_v,au_nw_v,au_se_v,au_ne_v &
 ,aw_sb_v,aw_nb_v,aw_st_v,aw_nt_v &
 , aw_w_w, aw_e_w, aw_s_w, aw_n_w &
 , aw_b_w, aw_t_w, aw_p_w         &
 ,au_bw_w,au_tw_w,au_be_w,au_te_w &
 ,av_bs_w,av_ts_w,av_bn_w,av_tn_w)

call cpy(ni,nj,nk,u,un)
call cpy(ni,nj,nk,v,vn)
call cpy(ni,nj,nk,w,wn)

call solu_sor4(ipara,ID,nID,ndiv,ni,nj,nk,key,sendjb,recvjb &
  ,nstep,imon_t,dt,rho,rhon &
  , au_w_u, au_e_u, au_s_u, au_n_u &
  , au_b_u, au_t_u, au_p_u         &
  ,av_ws_u,av_es_u,av_wn_u,av_en_u &
  ,aw_wb_u,aw_eb_u,aw_wt_u,aw_et_u &
  , av_w_v, av_e_v, av_s_v, av_n_v &
  , av_b_v, av_t_v, av_p_v         &
  ,au_sw_v,au_nw_v,au_se_v,au_ne_v &
  ,aw_sb_v,aw_nb_v,aw_st_v,aw_nt_v &
  , aw_w_w, aw_e_w, aw_s_w, aw_n_w &
  , aw_b_w, aw_t_w, aw_p_w         &
  ,au_bw_w,au_tw_w,au_be_w,au_te_w &
  ,av_bs_w,av_ts_w,av_bn_w,av_tn_w &
  ,src_u,src_v,src_w,un,vn,wn,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))

call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,un)
call bnd_periodic(ni,nj,nk,vn)
call bnd_periodic(ni,nj,nk,wn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,un)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wn)

call caldiv(ipara,ndiv,ni,nj,nk,dxinv,dyinv,dzinv,dt,un,vn,wn,div,err_div0)

call solp_fft_tdma1(ni,nj,nk,fftdata,div,phir_r,phir_i)
call trans_r2w(ipara,ID,ndiv,svall,key,sendbuf,recvbuf,phir_r,phiw_r)
call trans_r2w(ipara,ID,ndiv,svall,key,sendbuf,recvbuf,phir_i,phiw_i)
call solp_fft_tdma2(ID,ndiv,svall,rhog,dxinv,dzinv,as_p_fft,an_p_fft,ap_p_fft,atdma,btdma_r,btdma_i,phiw_r,phiw_i)
call trans_w2r(ipara,ID,ndiv,svall,key,sendbuf,recvbuf,phiw_r,phir_r)
call trans_w2r(ipara,ID,ndiv,svall,key,sendbuf,recvbuf,phiw_i,phir_i)
call solp_fft_tdma3(ni,nj,nk,fftdata,phir_r,phir_i,dp)
call bnd_avzero(ipara,ndiv,ni,nj,nk,dp)
call bnd_neumann(nID,ni,nj,nk,dp)
call bnd_periodic(ni,nj,nk,dp)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,dp)
call solp_fft_tdma4(ipara,ID,ndiv,ni,nj,nk,nstep,imon_t,rhog,dxinv,dyinv,dzinv,div,dp)

call corunp_explicit(nID,ni,nj,nk,rhog,dxinv,dyinv,dzinv,dt,dp,phat,un,vn,wn,pn)

call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1_array,l2_array,phi(-2,-2,-2,l))
call bnd_periodic(ni,nj,nk,un)
call bnd_periodic(ni,nj,nk,vn)
call bnd_periodic(ni,nj,nk,wn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,un)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wn)

call bnd_avzero(ipara,ndiv,ni,nj,nk,pn)
call bnd_neumann(nID,ni,nj,nk,pn)
call bnd_periodic(ni,nj,nk,pn)
call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,pn)

do l=1,nbub
call cal_center(ipara,ID,svall(2),ni,nj,nk,xl,yl,zl,dt,rhog,phi(-2,-2,-2,l) &
        ,center,center_pre1,center_pre2,area_pre1,area_pre2,lap,velocity)
call cal_vel(ipara,ni,nj,nk,particle_radius,rhol,rmul,u,v,w,phi(-2,-2,-2,l),vel,re)
if(mod(nstep,imkuvp).eq.0)then
if(ID.eq.0)then
call mk_center(nstep,time,center,velocity,l)
call mk_vel(nstep,time,vel,re,velocity,l)
endif
endif
enddo

if(mod(nstep,imon_t).eq.0)then
  call caldiv(ipara,ndiv,ni,nj,nk,dxinv,dyinv,dzinv,dt,un,vn,wn,div,err_div)
  if(ID.eq.0)then
    write(*,'("Err_div  ",1i10,20e20.10)') &
      nstep &
     ,err_div/err_div0 &
     ,err_div &
     ,err_div0
   endif
endif

call cpy(ni,nj,nk,u   ,uo  )
call cpy(ni,nj,nk,v   ,vo  )
call cpy(ni,nj,nk,w   ,wo  )
call cpy(ni,nj,nk,p   ,po  )
call cpy(ni,nj,nk,un  ,u   )
call cpy(ni,nj,nk,vn  ,v   )
call cpy(ni,nj,nk,wn  ,w   )
call cpy(ni,nj,nk,pn  ,p   )

do l=0,nbub
call cpy(ni,nj,nk,phin(-2,-2,-2,l),phi(-2,-2,-2,l))
enddo

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)

if(mod(nstep,idout).eq.0)then                                               
  call dataou(ipara,ID,ni,nj,nk,nbub,nstep,time,u,v,w,p,uo,vo,wo,po,phi)
endif

if(mod(nstep,imkuvp).eq.0)then
  j = particle_radius/dy
  call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,0,1,1,wrk_all,  u,  u_slice)
  call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,1,1,0,wrk_all,  w,  w_slice)
  call mk_slice_j(ipara,ID,nID,ndiv,ni,nj,nk,j,1,1,1,wrk_all,phi,phi_slice)
  if(ID.eq.0)then
    include'mkuwp.h'
  endif
endif


if(mod(nstep,imkvtk).eq.0)then
  call calvorx(ni,nj,nk,dyinv,dzinv,v,w,vorx)
  call bnd_periodic(ni,nj,nk,vorx)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vorx)
  call cal_advu(ni,nj,nk,dxinv,dyinv,dzinv,u,v,w,adv_u,adv_v,adv_w)
  call calq(ni,nj,nk,dxinv,dyinv,dzinv,adv_u,adv_v,adv_w,q)
  call bnd_periodic(ni,nj,nk,q)
  call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,q)
  do l=1,nbub
  call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all, phi(-2,-2,-2,l), phil_all(-2,-2,-2,l)) 
  enddo
  call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all, phi, phi_all)
  call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,   q,   q_all)
  call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,vorx,vorx_all)
  call mk_all(ipara,ID,nID,ndiv,ni,nj,nk,wrk_all,   p,   p_all)
  if(ID.eq.0)then
  !do l=1,nbub
  !call mkvtk_phil(svall,nstep,dx,dy,dz,phil_all(-2,-2,-2,l),l)
  !enddo

  call mkvtk_phi(svall,nstep,dx,dy,dz, phi_all)
  call mkvtk_p(svall,nstep,dx,dy,dz,   p_all)
  call find_interface_positions(ni, nj, nk, phi_all, dx, dy, xl)
  !  call   mkvtk_q(svall,nstep,dx,dy,dz,vorx_all,q_all)
  endif
endif

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)

enddo !nstep

call mpi_barrier(mpi_comm_world,ierr)
call flush(6)
call mpi_finalize(ierr)

stop
end program main
