      subroutine solu_sor4(ipara,ID,nID,ndiv,ni,nj,nk,key,sendjb,recvjb
     & ,nstep,imon_t
     & ,dt
     & ,rho,rhon
     & , au_w_u, au_e_u, au_s_u, au_n_u
     & , au_b_u, au_t_u, au_p_u
     & ,av_ws_u,av_es_u,av_wn_u,av_en_u
     & ,aw_wb_u,aw_eb_u,aw_wt_u,aw_et_u
     & , av_w_v, av_e_v, av_s_v, av_n_v
     & , av_b_v, av_t_v, av_p_v
     & ,au_sw_v,au_nw_v,au_se_v,au_ne_v
     & ,aw_sb_v,aw_nb_v,aw_st_v,aw_nt_v
     & , aw_w_w, aw_e_w, aw_s_w, aw_n_w
     & , aw_b_w, aw_t_w, aw_p_w
     & ,au_bw_w,au_tw_w,au_be_w,au_te_w
     & ,av_bs_w,av_ts_w,av_bn_w,av_tn_w
     & ,src_u,src_v,src_w,un,vn,wn,uwall,dy,l1)

      implicit none
      include 'mpif.h'
      include 'param.h'
      integer ipara,ID,nID(6),ndiv
      integer ni,nj,nk
      integer  key(2,2)
      real*8 recvjb(-2:ni+3,-2:nk+3,3,2)
      real*8 sendjb(-2:ni+3,-2:nk+3,3,2)
      integer nstep,imon_t
      real*8 dt
      real*8     rho(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8    rhon(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_w_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_e_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_s_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_n_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_b_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_t_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_p_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_ws_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_es_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_wn_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_en_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_wb_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_eb_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_wt_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_et_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_w_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_e_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_s_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_n_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_b_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_t_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_p_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_sw_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_nw_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_se_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_ne_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_sb_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_nb_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_st_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_nt_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_w_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_e_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_s_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_n_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_b_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_t_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_p_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_bw_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_tw_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_be_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_te_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_bs_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_ts_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_bn_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_tn_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_u(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_v(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8   src_w(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      un(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      vn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8      wn(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 uwall,dy,l1

      integer ierr
      real*8 sendbuf(1),recvbuf(1)
      integer i,j,k,jen,lp,lsor4
      integer isor(8),jsor(8),ksor(8)
      real*8 dtinv,rho00
      real*8 err,err0,du,dv,dw,visn

      isor(1)=1
      jsor(1)=1
      ksor(1)=1
      isor(2)=2
      jsor(2)=2
      ksor(2)=2

      isor(3)=1
      jsor(3)=1
      ksor(3)=2
      isor(4)=2
      jsor(4)=2
      ksor(4)=1

      isor(5)=1
      jsor(5)=2
      ksor(5)=1
      isor(6)=2
      jsor(6)=1
      ksor(6)=2

      isor(7)=2
      jsor(7)=1
      ksor(7)=1
      isor(8)=1
      jsor(8)=2
      ksor(8)=2

      dtinv=1.0d0/dt

      jen=nj
      if(nID(Y_PLUS).lt.0)jen=nj-1

cc
cc<sor4
cc
      do lp=0,20
      err=0.0d0

      do lsor4=1,8
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(visn,rho00,du)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(isor,jsor,ksor,lsor4)
!$OMP$ SHARED(rho,rhon)
!$OMP$ SHARED( au_w_u, au_e_u, au_s_u, au_n_u)
!$OMP$ SHARED( au_b_u, au_t_u, au_p_u)
!$OMP$ SHARED(av_ws_u,av_es_u,av_wn_u,av_en_u)
!$OMP$ SHARED(aw_wb_u,aw_eb_u,aw_wt_u,aw_et_u)
!$OMP$ SHARED(un,vn,wn,src_u,dtinv)
!$OMP$ REDUCTION(+:err)
      do k=ksor(lsor4),nk,2
      do j=jsor(lsor4),nj,2
      do i=isor(lsor4),ni,2
      rho00=(
     1   rho(i,j,k)+ rho(i+1,j,k)
     2 +rhon(i,j,k)+rhon(i+1,j,k))*0.25d0
      visn=(
     & + au_b_u(i,j,k)*un(i  ,j  ,k-1)
     & + au_s_u(i,j,k)*un(i  ,j-1,k  )
     & + au_w_u(i,j,k)*un(i-1,j  ,k  )
     & + au_e_u(i,j,k)*un(i+1,j  ,k  )
     & + au_n_u(i,j,k)*un(i  ,j+1,k  )
     & + au_t_u(i,j,k)*un(i  ,j  ,k+1)
     & +av_ws_u(i,j,k)*vn(i  ,j-1,k  )
     & +av_es_u(i,j,k)*vn(i+1,j-1,k  )
     & +av_wn_u(i,j,k)*vn(i  ,j  ,k  )
     & +av_en_u(i,j,k)*vn(i+1,j  ,k  )
     & +aw_wb_u(i,j,k)*wn(i  ,j  ,k-1)
     & +aw_eb_u(i,j,k)*wn(i+1,j  ,k-1)
     & +aw_wt_u(i,j,k)*wn(i  ,j  ,k  )
     & +aw_et_u(i,j,k)*wn(i+1,j  ,k  ) )/rho00
      du=(visn*0.5d0+src_u(i,j,k))/(dtinv+au_p_u(i,j,k)/rho00*0.5d0)
     & -un(i,j,k)
      un(i,j,k)=un(i,j,k)+du
      err=err+du**2
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      if(mod(lsor4,2).eq.0)then
      call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1)
      call bnd_periodic(ni,nj,nk,un)
      call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,un)
      endif
      enddo

      do lsor4=1,8
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(visn,rho00,dv)
!$OMP$ SHARED(ni,jen,nk,ID,nj)
!$OMP$ SHARED(isor,jsor,ksor,lsor4)
!$OMP$ SHARED(rho,rhon)
!$OMP$ SHARED( av_w_v, av_e_v, av_s_v, av_n_v)
!$OMP$ SHARED( av_b_v, av_t_v, av_p_v)
!$OMP$ SHARED(au_sw_v,au_nw_v,au_se_v,au_ne_v)
!$OMP$ SHARED(aw_sb_v,aw_nb_v,aw_st_v,aw_nt_v)
!$OMP$ SHARED(un,vn,wn,src_v,dtinv)
!$OMP$ REDUCTION(+:err)
      do k=ksor(lsor4),nk ,2
      do j=jsor(lsor4),jen,2
      do i=isor(lsor4),ni ,2
      rho00=(
     1   rho(i,j,k)+ rho(i,j+1,k)
     2 +rhon(i,j,k)+rhon(i,j+1,k))*0.25d0
      visn=(
     & + av_b_v(i,j,k)*vn(i  ,j  ,k-1)
     & + av_s_v(i,j,k)*vn(i  ,j-1,k  )
     & + av_w_v(i,j,k)*vn(i-1,j  ,k  )
     & + av_e_v(i,j,k)*vn(i+1,j  ,k  )
     & + av_n_v(i,j,k)*vn(i  ,j+1,k  )
     & + av_t_v(i,j,k)*vn(i  ,j  ,k+1)
     & +au_sw_v(i,j,k)*un(i-1,j  ,k  )
     & +au_se_v(i,j,k)*un(i  ,j  ,k  )
     & +au_nw_v(i,j,k)*un(i-1,j+1,k  )
     & +au_ne_v(i,j,k)*un(i  ,j+1,k  )
     & +aw_sb_v(i,j,k)*wn(i  ,j  ,k-1)
     & +aw_nb_v(i,j,k)*wn(i  ,j+1,k-1)
     & +aw_st_v(i,j,k)*wn(i  ,j  ,k  )
     & +aw_nt_v(i,j,k)*wn(i  ,j+1,k  ) )/rho00
      dv=(visn*0.5d0+src_v(i,j,k))/(dtinv+av_p_v(i,j,k)/rho00*0.5d0)
     & -vn(i,j,k)
      vn(i,j,k)=vn(i,j,k)+dv
      err=err+dv**2
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      if(mod(lsor4,2).eq.0)then
      call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1)
      call bnd_periodic(ni,nj,nk,vn)
      call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,vn)
      endif
      enddo

      do lsor4=1,8
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(visn,rho00,dw)
!$OMP$ SHARED(ni,nj,nk,ID)
!$OMP$ SHARED(isor,jsor,ksor,lsor4)
!$OMP$ SHARED(rho,rhon)
!$OMP$ SHARED( aw_w_w, aw_e_w, aw_s_w, aw_n_w)
!$OMP$ SHARED( aw_b_w, aw_t_w, aw_p_w)
!$OMP$ SHARED(au_bw_w,au_tw_w,au_be_w,au_te_w)
!$OMP$ SHARED(av_bs_w,av_ts_w,av_bn_w,av_tn_w)
!$OMP$ SHARED(un,vn,wn,src_w,dtinv)
!$OMP$ REDUCTION(+:err)
      do k=ksor(lsor4),nk,2
      do j=jsor(lsor4),nj,2
      do i=isor(lsor4),ni,2
      rho00=(
     1   rho(i,j,k)+ rho(i,j,k+1)
     2 +rhon(i,j,k)+rhon(i,j,k+1))*0.25d0
      visn=(
     & + aw_b_w(i,j,k)*wn(i  ,j  ,k-1)
     & + aw_s_w(i,j,k)*wn(i  ,j-1,k  )
     & + aw_w_w(i,j,k)*wn(i-1,j  ,k  )
     & + aw_e_w(i,j,k)*wn(i+1,j  ,k  )
     & + aw_n_w(i,j,k)*wn(i  ,j+1,k  )
     & + aw_t_w(i,j,k)*wn(i  ,j  ,k+1)
     & +au_bw_w(i,j,k)*un(i-1,j  ,k  )
     & +au_be_w(i,j,k)*un(i  ,j  ,k  )
     & +au_tw_w(i,j,k)*un(i-1,j  ,k+1)
     & +au_te_w(i,j,k)*un(i  ,j  ,k+1)
     & +av_bs_w(i,j,k)*vn(i  ,j-1,k  )
     & +av_bn_w(i,j,k)*vn(i  ,j  ,k  )
     & +av_ts_w(i,j,k)*vn(i  ,j-1,k+1)
     & +av_tn_w(i,j,k)*vn(i  ,j  ,k+1) )/rho00
      dw=(visn*0.5d0+src_w(i,j,k))/(dtinv+aw_p_w(i,j,k)/rho00*0.5d0)
     & -wn(i,j,k)
      wn(i,j,k)=wn(i,j,k)+dw
      err=err+dw**2
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO
      if(mod(lsor4,2).eq.0)then
      call bndu(nID,ni,nj,nk,un,vn,wn,uwall,dy,l1)
      call bnd_periodic(ni,nj,nk,wn)
      call bnd_comm(ipara,nID,ni,nj,nk,key,sendjb,recvjb,wn)
      endif
      enddo

      if(ipara.eq.1)then
      sendbuf(1) = err
      recvbuf(1) = err
      call mpi_allreduce(sendbuf,recvbuf,1,MPI_REAL8,MPI_SUM
     & ,MPI_COMM_WORLD,ierr)
      err = recvbuf(1)
      endif

      err=sqrt(err/dble(ni*nj*nk*ndiv))+1.0d-99
      if(lp.eq.0)err0=err
      if(err/err0.le.1.0d-3)goto 100
      enddo

 100  continue

      if(mod(nstep,imon_t).eq.0.and.ID.eq.0)then
      write(*,'("Err_solu ",2i10,20e20.10)')
     2  nstep
     3 ,lp
     4 ,err/err0
     5 ,err
     6 ,err0
      endif

      return
      end
