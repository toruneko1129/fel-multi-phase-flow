ccc
ccc<compute expansion coefficients for determining the viscous terms
ccc from the viscosity rmuk and the density rhok
ccc 

      subroutine cal_arith_coef_vis(ni,nj,nk,dxinv,dyinv,dzinv
     & ,    rmuk
     & , au_w_uk, au_e_uk, au_s_uk, au_n_uk
     & , au_b_uk, au_t_uk, au_p_uk
     & ,av_ws_uk,av_es_uk,av_wn_uk,av_en_uk
     & ,aw_wb_uk,aw_eb_uk,aw_wt_uk,aw_et_uk
     & , av_w_vk, av_e_vk, av_s_vk, av_n_vk
     & , av_b_vk, av_t_vk, av_p_vk
     & ,au_sw_vk,au_nw_vk,au_se_vk,au_ne_vk
     & ,aw_sb_vk,aw_nb_vk,aw_st_vk,aw_nt_vk
     & , aw_w_wk, aw_e_wk, aw_s_wk, aw_n_wk
     & , aw_b_wk, aw_t_wk, aw_p_wk
     & ,au_bw_wk,au_tw_wk,au_be_wk,au_te_wk
     & ,av_bs_wk,av_ts_wk,av_bn_wk,av_tn_wk)

      implicit none
      integer ni,nj,nk
      real*8 dxinv,dyinv,dzinv
      real*8     rmuk(-2:ni+3,-2:nj+3,-2:nk+3)

      real*8  au_w_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_e_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_s_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_n_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_b_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_t_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  au_p_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_ws_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_es_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_wn_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_en_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_wb_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_eb_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_wt_uk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_et_uk(-2:ni+3,-2:nj+3,-2:nk+3)

      real*8  av_w_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_e_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_s_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_n_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_b_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_t_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  av_p_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_sw_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_nw_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_se_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_ne_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_sb_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_nb_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_st_vk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 aw_nt_vk(-2:ni+3,-2:nj+3,-2:nk+3)

      real*8  aw_w_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_e_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_s_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_n_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_b_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_t_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8  aw_p_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_bw_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_tw_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_be_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 au_te_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_bs_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_ts_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_bn_wk(-2:ni+3,-2:nj+3,-2:nk+3)
      real*8 av_tn_wk(-2:ni+3,-2:nj+3,-2:nk+3)

      integer i,j,k
      real*8 rmu_xm,rmu_ym,rmu_zm,rmu_xp,rmu_yp,rmu_zp


ccc
ccc<initially clear the coefficients
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED( au_w_uk, au_e_uk, au_s_uk, au_n_uk)
!$OMP$ SHARED( au_b_uk, au_t_uk, au_p_uk)
!$OMP$ SHARED(av_ws_uk,av_es_uk,av_wn_uk,av_en_uk)
!$OMP$ SHARED(aw_wb_uk,aw_eb_uk,aw_wt_uk,aw_et_uk)
!$OMP$ SHARED( av_w_vk, av_e_vk, av_s_vk, av_n_vk)
!$OMP$ SHARED( av_b_vk, av_t_vk, av_p_vk)
!$OMP$ SHARED(au_sw_vk,au_nw_vk,au_se_vk,au_ne_vk)
!$OMP$ SHARED(aw_sb_vk,aw_nb_vk,aw_st_vk,aw_nt_vk)
!$OMP$ SHARED( aw_w_wk, aw_e_wk, aw_s_wk, aw_n_wk)
!$OMP$ SHARED( aw_b_wk, aw_t_wk, aw_p_wk)
!$OMP$ SHARED(au_bw_wk,au_tw_wk,au_be_wk,au_te_wk)
!$OMP$ SHARED(av_bs_wk,av_ts_wk,av_bn_wk,av_tn_wk)
      do k=-2,nk+3
      do j=-2,nj+3
      do i=-2,ni+3
       au_w_uk(i,j,k)=0.0d0
       au_e_uk(i,j,k)=0.0d0
       au_s_uk(i,j,k)=0.0d0
       au_n_uk(i,j,k)=0.0d0
       au_b_uk(i,j,k)=0.0d0
       au_t_uk(i,j,k)=0.0d0
       au_p_uk(i,j,k)=0.0d0
      av_ws_uk(i,j,k)=0.0d0
      av_es_uk(i,j,k)=0.0d0
      av_wn_uk(i,j,k)=0.0d0
      av_en_uk(i,j,k)=0.0d0
      aw_wb_uk(i,j,k)=0.0d0
      aw_eb_uk(i,j,k)=0.0d0
      aw_wt_uk(i,j,k)=0.0d0
      aw_et_uk(i,j,k)=0.0d0
       av_w_vk(i,j,k)=0.0d0
       av_e_vk(i,j,k)=0.0d0
       av_s_vk(i,j,k)=0.0d0
       av_n_vk(i,j,k)=0.0d0
       av_b_vk(i,j,k)=0.0d0
       av_t_vk(i,j,k)=0.0d0
       av_p_vk(i,j,k)=0.0d0
      au_sw_vk(i,j,k)=0.0d0
      au_nw_vk(i,j,k)=0.0d0
      au_se_vk(i,j,k)=0.0d0
      au_ne_vk(i,j,k)=0.0d0
      aw_sb_vk(i,j,k)=0.0d0
      aw_nb_vk(i,j,k)=0.0d0
      aw_st_vk(i,j,k)=0.0d0
      aw_nt_vk(i,j,k)=0.0d0
       aw_w_wk(i,j,k)=0.0d0
       aw_e_wk(i,j,k)=0.0d0
       aw_s_wk(i,j,k)=0.0d0
       aw_n_wk(i,j,k)=0.0d0
       aw_b_wk(i,j,k)=0.0d0
       aw_t_wk(i,j,k)=0.0d0
       aw_p_wk(i,j,k)=0.0d0
      au_bw_wk(i,j,k)=0.0d0
      au_tw_wk(i,j,k)=0.0d0
      au_be_wk(i,j,k)=0.0d0
      au_te_wk(i,j,k)=0.0d0
      av_bs_wk(i,j,k)=0.0d0
      av_ts_wk(i,j,k)=0.0d0
      av_bn_wk(i,j,k)=0.0d0
      av_tn_wk(i,j,k)=0.0d0
      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<set coefficients for the viscous term vis_uk
ccc       vis_uk(i,j,k)=(
ccc  & + au_b_uk(i,j,k)*uk(i  ,j  ,k-1)
ccc  & + au_s_uk(i,j,k)*uk(i  ,j-1,k  )
ccc  & + au_w_uk(i,j,k)*uk(i-1,j  ,k  )
ccc  & - au_p_uk(i,j,k)*uk(i  ,j  ,k  )
ccc  & + au_e_uk(i,j,k)*uk(i+1,j  ,k  )
ccc  & + au_n_uk(i,j,k)*uk(i  ,j+1,k  )
ccc  & + au_t_uk(i,j,k)*uk(i  ,j  ,k+1)
ccc  & +av_ws_uk(i,j,k)*vk(i  ,j-1,k  )
ccc  & +av_es_uk(i,j,k)*vk(i+1,j-1,k  )
ccc  & +av_wn_uk(i,j,k)*vk(i  ,j  ,k  )
ccc  & +av_en_uk(i,j,k)*vk(i+1,j  ,k  )
ccc  & +aw_wb_uk(i,j,k)*wk(i  ,j  ,k-1)
ccc  & +aw_eb_uk(i,j,k)*wk(i+1,j  ,k-1)
ccc  & +aw_wt_uk(i,j,k)*wk(i  ,j  ,k  )
ccc  & +aw_et_uk(i,j,k)*wk(i+1,j  ,k  )
ccc  & )
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(rmu_xm,rmu_ym,rmu_zm,rmu_xp,rmu_yp,rmu_zp)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rmuk,dxinv,dyinv,dzinv)
!$OMP$ SHARED( au_w_uk, au_e_uk, au_s_uk, au_n_uk)
!$OMP$ SHARED( au_b_uk, au_t_uk, au_p_uk)
!$OMP$ SHARED(av_ws_uk,av_es_uk,av_wn_uk,av_en_uk)
!$OMP$ SHARED(aw_wb_uk,aw_eb_uk,aw_wt_uk,aw_et_uk)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      rmu_xm= rmuk(i  ,j  ,k  )
      rmu_xp= rmuk(i+1,j  ,k  )
      rmu_ym=0.25d0*(
     1       +rmuk(i  ,j-1,k  )
     2       +rmuk(i+1,j-1,k  )
     3       +rmuk(i  ,j  ,k  )
     4       +rmuk(i+1,j  ,k  ))
      rmu_yp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i+1,j  ,k  )
     3       +rmuk(i  ,j+1,k  )
     4       +rmuk(i+1,j+1,k  ))
      rmu_zm=0.25d0*(
     1       +rmuk(i  ,j  ,k-1)
     2       +rmuk(i+1,j  ,k-1)
     3       +rmuk(i  ,j  ,k  )
     4       +rmuk(i+1,j  ,k  ))
      rmu_zp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i+1,j  ,k  )
     3       +rmuk(i  ,j  ,k+1)
     4       +rmuk(i+1,j  ,k+1))

        au_b_uk(i,j,k)= dzinv**2*rmu_zm
        au_s_uk(i,j,k)= dyinv**2*rmu_ym
        au_w_uk(i,j,k)= dxinv**2*rmu_xm*2.0d0
        au_e_uk(i,j,k)= dxinv**2*rmu_xp*2.0d0
        au_n_uk(i,j,k)= dyinv**2*rmu_yp
        au_t_uk(i,j,k)= dzinv**2*rmu_zp
       av_ws_uk(i,j,k)=+dxinv*dyinv*rmu_ym
       av_es_uk(i,j,k)=-dxinv*dyinv*rmu_ym
       av_wn_uk(i,j,k)=-dxinv*dyinv*rmu_yp
       av_en_uk(i,j,k)=+dxinv*dyinv*rmu_yp
       aw_wb_uk(i,j,k)=+dxinv*dzinv*rmu_zm
       aw_eb_uk(i,j,k)=-dxinv*dzinv*rmu_zm
       aw_wt_uk(i,j,k)=-dxinv*dzinv*rmu_zp
       aw_et_uk(i,j,k)=+dxinv*dzinv*rmu_zp
        au_p_uk(i,j,k)=
     & +au_b_uk(i,j,k)
     & +au_s_uk(i,j,k)
     & +au_w_uk(i,j,k)
     & +au_e_uk(i,j,k)
     & +au_n_uk(i,j,k)
     & +au_t_uk(i,j,k)

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<set coefficients for the viscous term vis_vk
ccc       vis_vk(i,j,k)=(
ccc  & + av_b_vk(i,j,k)*vk(i  ,j  ,k-1)
ccc  & + av_s_vk(i,j,k)*vk(i  ,j-1,k  )
ccc  & + av_w_vk(i,j,k)*vk(i-1,j  ,k  )
ccc  & - av_p_vk(i,j,k)*vk(i  ,j  ,k  )
ccc  & + av_e_vk(i,j,k)*vk(i+1,j  ,k  )
ccc  & + av_n_vk(i,j,k)*vk(i  ,j+1,k  )
ccc  & + av_t_vk(i,j,k)*vk(i  ,j  ,k+1)
ccc  & +au_sw_vk(i,j,k)*uk(i-1,j  ,k  )
ccc  & +au_se_vk(i,j,k)*uk(i  ,j  ,k  )
ccc  & +au_nw_vk(i,j,k)*uk(i-1,j+1,k  )
ccc  & +au_ne_vk(i,j,k)*uk(i  ,j+1,k  )
ccc  & +aw_sb_vk(i,j,k)*wk(i  ,j  ,k-1)
ccc  & +aw_nb_vk(i,j,k)*wk(i  ,j+1,k-1)
ccc  & +aw_st_vk(i,j,k)*wk(i  ,j  ,k  )
ccc  & +aw_nt_vk(i,j,k)*wk(i  ,j+1,k  )
ccc  & )
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(rmu_xm,rmu_ym,rmu_zm,rmu_xp,rmu_yp,rmu_zp)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rmuk,dxinv,dyinv,dzinv)
!$OMP$ SHARED( av_w_vk, av_e_vk, av_s_vk, av_n_vk)
!$OMP$ SHARED( av_b_vk, av_t_vk, av_p_vk)
!$OMP$ SHARED(au_sw_vk,au_nw_vk,au_se_vk,au_ne_vk)
!$OMP$ SHARED(aw_sb_vk,aw_nb_vk,aw_st_vk,aw_nt_vk)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      rmu_xm=0.25d0*(
     1       +rmuk(i-1,j  ,k  )
     2       +rmuk(i  ,j  ,k  )
     3       +rmuk(i-1,j+1,k  )
     4       +rmuk(i  ,j+1,k  ))
      rmu_xp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i+1,j  ,k  )
     3       +rmuk(i  ,j+1,k  )
     4       +rmuk(i+1,j+1,k  ))
      rmu_ym= rmuk(i  ,j  ,k  )
      rmu_yp= rmuk(i  ,j+1,k  )
      rmu_zm=0.25d0*(
     1       +rmuk(i  ,j  ,k-1)
     2       +rmuk(i  ,j+1,k-1)
     3       +rmuk(i  ,j  ,k  )
     4       +rmuk(i  ,j+1,k  ))
      rmu_zp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i  ,j+1,k  )
     3       +rmuk(i  ,j  ,k+1)
     4       +rmuk(i  ,j+1,k+1))

        av_b_vk(i,j,k)= dzinv**2*rmu_zm
        av_s_vk(i,j,k)= dyinv**2*rmu_ym*2.0d0
        av_w_vk(i,j,k)= dxinv**2*rmu_xm
        av_e_vk(i,j,k)= dxinv**2*rmu_xp
        av_n_vk(i,j,k)= dyinv**2*rmu_yp*2.0d0
        av_t_vk(i,j,k)= dzinv**2*rmu_zp
       au_sw_vk(i,j,k)=+dxinv*dyinv*rmu_xm
       au_nw_vk(i,j,k)=-dxinv*dyinv*rmu_xm
       au_se_vk(i,j,k)=-dxinv*dyinv*rmu_xp
       au_ne_vk(i,j,k)=+dxinv*dyinv*rmu_xp
       aw_sb_vk(i,j,k)=+dyinv*dzinv*rmu_zm
       aw_nb_vk(i,j,k)=-dyinv*dzinv*rmu_zm
       aw_st_vk(i,j,k)=-dyinv*dzinv*rmu_zp
       aw_nt_vk(i,j,k)=+dyinv*dzinv*rmu_zp
        av_p_vk(i,j,k)=
     & +av_b_vk(i,j,k)
     & +av_s_vk(i,j,k)
     & +av_w_vk(i,j,k)
     & +av_e_vk(i,j,k)
     & +av_n_vk(i,j,k)
     & +av_t_vk(i,j,k)

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

ccc
ccc<set coefficients for the viscous term vis_wk
ccc       vis_wk(i,j,k)=(
ccc  & + aw_b_wk(i,j,k)*wk(i  ,j  ,k-1)
ccc  & + aw_s_wk(i,j,k)*wk(i  ,j-1,k  )
ccc  & + aw_w_wk(i,j,k)*wk(i-1,j  ,k  )
ccc  & - aw_p_wk(i,j,k)*wk(i  ,j  ,k  )
ccc  & + aw_e_wk(i,j,k)*wk(i+1,j  ,k  )
ccc  & + aw_n_wk(i,j,k)*wk(i  ,j+1,k  )
ccc  & + aw_t_wk(i,j,k)*wk(i  ,j  ,k+1)
ccc  & +au_bw_wk(i,j,k)*uk(i-1,j  ,k  )
ccc  & +au_be_wk(i,j,k)*uk(i  ,j  ,k  )
ccc  & +au_tw_wk(i,j,k)*uk(i-1,j  ,k+1)
ccc  & +au_te_wk(i,j,k)*uk(i  ,j  ,k+1)
ccc  & +av_bs_wk(i,j,k)*vk(i  ,j-1,k  )
ccc  & +av_bn_wk(i,j,k)*vk(i  ,j  ,k  )
ccc  & +av_ts_wk(i,j,k)*vk(i  ,j-1,k+1)
ccc  & +av_tn_wk(i,j,k)*vk(i  ,j  ,k+1)
ccc  & )
ccc

!$OMP  PARALLEL DO
!$OMP$ SCHEDULE(static,1)
!$OMP$ DEFAULT(none)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ PRIVATE(rmu_xm,rmu_ym,rmu_zm,rmu_xp,rmu_yp,rmu_zp)
!$OMP$ SHARED(ni,nj,nk)
!$OMP$ SHARED(rmuk,dxinv,dyinv,dzinv)
!$OMP$ SHARED( aw_w_wk, aw_e_wk, aw_s_wk, aw_n_wk)
!$OMP$ SHARED( aw_b_wk, aw_t_wk, aw_p_wk)
!$OMP$ SHARED(au_bw_wk,au_tw_wk,au_be_wk,au_te_wk)
!$OMP$ SHARED(av_bs_wk,av_ts_wk,av_bn_wk,av_tn_wk)
      do k=-1,nk+2
      do j=-1,nj+2
      do i=-1,ni+2

      rmu_xm=0.25d0*(
     1       +rmuk(i-1,j  ,k  )
     2       +rmuk(i  ,j  ,k  )
     3       +rmuk(i-1,j  ,k+1)
     4       +rmuk(i  ,j  ,k+1))
      rmu_xp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i+1,j  ,k  )
     3       +rmuk(i  ,j  ,k+1)
     4       +rmuk(i+1,j  ,k+1))
      rmu_ym=0.25d0*(
     1       +rmuk(i  ,j-1,k  )
     2       +rmuk(i  ,j  ,k  )
     3       +rmuk(i  ,j-1,k+1)
     4       +rmuk(i  ,j  ,k+1))
      rmu_yp=0.25d0*(
     1       +rmuk(i  ,j  ,k  )
     2       +rmuk(i  ,j+1,k  )
     3       +rmuk(i  ,j  ,k+1)
     4       +rmuk(i  ,j+1,k+1))
      rmu_zm= rmuk(i  ,j  ,k  )
      rmu_zp= rmuk(i  ,j  ,k+1)

        aw_b_wk(i,j,k)= dzinv**2*rmu_zm*2.0d0
        aw_s_wk(i,j,k)= dyinv**2*rmu_ym
        aw_w_wk(i,j,k)= dxinv**2*rmu_xm
        aw_e_wk(i,j,k)= dxinv**2*rmu_xp
        aw_n_wk(i,j,k)= dyinv**2*rmu_yp
        aw_t_wk(i,j,k)= dzinv**2*rmu_zp*2.0d0
       au_bw_wk(i,j,k)=+dxinv*dzinv*rmu_xm
       au_tw_wk(i,j,k)=-dxinv*dzinv*rmu_xm
       au_be_wk(i,j,k)=-dxinv*dzinv*rmu_xp
       au_te_wk(i,j,k)=+dxinv*dzinv*rmu_xp
       av_bs_wk(i,j,k)=+dyinv*dzinv*rmu_ym
       av_ts_wk(i,j,k)=-dyinv*dzinv*rmu_ym
       av_bn_wk(i,j,k)=-dyinv*dzinv*rmu_yp
       av_tn_wk(i,j,k)=+dyinv*dzinv*rmu_yp
        aw_p_wk(i,j,k)=
     & +aw_b_wk(i,j,k)
     & +aw_s_wk(i,j,k)
     & +aw_w_wk(i,j,k)
     & +aw_e_wk(i,j,k)
     & +aw_n_wk(i,j,k)
     & +aw_t_wk(i,j,k)

      enddo
      enddo
      enddo
!$OMP  END PARALLEL DO

      return
      end
