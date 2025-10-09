!<impose boundary conditions on the velocity components uk,vk,wk
!<localised navier slip
!<l1, l2: slip length of fluid[1,2]

subroutine bndu_localized_free(nID,ni,nj,nk,uk,vk,wk,uwall,dx,dy,l1,l2,phi)
  implicit none
  include 'param.h'

  integer nID(6), ni,nj,nk
  real*8  uk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  vk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  wk(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  uwall, dx, dy, epsb
  real*8  l1(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  l2(-2:ni+3,-2:nj+3,-2:nk+3)
  real*8  phi(-2:ni+3,-2:nj+3,-2:nk+3)

  !==== 追加の型宣言（implicit none 対応）=========================
  integer :: i,k,j        ! j を明示宣言（後半の周期コピーで使用）
  real*8  bell_weight           ! ← 関数の戻り値型を明示
  external bell_weight          ! ← 明示インターフェースが無い場合は付けると安全
  !==============================================================

  real*8  :: phi_av, inv_ls1, inv_ls2, ls_mix, ls_eff, coef1, coef2
  real*8  :: x_i, d, w_up, w_dn, w, tiny
  real*8  :: xup_top(-2:nk+3), xdn_top(-2:nk+3), xup_bot(-2:nk+3), xdn_bot(-2:nk+3)
  logical :: hup_top(-2:nk+3), hdn_top(-2:nk+3), hup_bot(-2:nk+3), hdn_bot(-2:nk+3)
  real*8  :: l_free, r0, r1, r2
  real*8  :: slip_len_profile
  external   slip_len_profile

  tiny  = 1.0d-300
  epsb  = dx * 2.5d0
  l_free = 1.0d0*dy          ! free-slip 相当の大滑り長さ
  r0     = 2.0d0*dx            ! コア半径
  r1     = 5.0d0*dx            ! 遷移外縁1
  r2     = 64.0d0*dx           ! 遷移外縁2

  ! 各壁の2交差を取得
  call find_two_contacts_on_wall(ni,nj,nk,phi,dx,  1, xup_bot,hup_bot, xdn_bot,hdn_bot)
  call find_two_contacts_on_wall(ni,nj,nk,phi,dx, nj, xup_top,hup_top, xdn_top,hdn_top)

  !================ 下壁 j=0 側 =================
  if (nID(Y_MINUS) .lt. 0) then
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,k,phi_av,inv_ls1,inv_ls2,ls_mix,ls_eff,coef1,coef2,x_i,d,w_up,w_dn,w) &
!$OMP& SHARED(ni,nk,dx,dy,uk,vk,wk,phi,uwall,l1,l2,epsb,xup_bot,xdn_bot,hup_bot,hdn_bot,tiny,l_free,r0,r1,r2)
    do k=-2, nk+3
      do i=-2, ni+3
        phi_av = 0.5d0*(phi(i,1,k) + phi(i+1,1,k))
        inv_ls1 =  phi_av        / max(l1(i,1,k), tiny)
        inv_ls2 = (1.0d0-phi_av) / max(l2(i,1,k), tiny)
        ls_mix  = 1.0d0 / max(inv_ls1 + inv_ls2, tiny)

        x_i = (dble(i)+0.5d0)*dx
        w_up = 0.0d0; w_dn = 0.0d0
        if (hup_bot(k)) then
          d = x_i - xup_bot(k)
          w_up = slip_len_profile(d, r0, r1, r2, l_free, ls_mix)
        end if
        if (hdn_bot(k)) then
          d = x_i - xdn_bot(k)
          w_dn = slip_len_profile(d, r0, r1, r2, l_free, ls_mix)
        end if

        ls_eff = max(w_up, w_dn)

        coef1 = (2.0d0*dy)         / (2.0d0*ls_eff + dy)
        coef2 = (2.0d0*ls_eff - dy) / (2.0d0*ls_eff + dy)

        uk(i,   0,k) = coef1*(-uwall) + coef2*uk(i,   1,k)
        uk(i,  -1,k) = 2.0d0*uk(i,   0,k) - uk(i,   1,k)
        uk(i,  -2,k) = 2.0d0*uk(i,  -1,k) - uk(i,   0,k)

        vk(i,  -2,k) = -vk(i,   2,k)
        vk(i,  -1,k) = -vk(i,   1,k)
        vk(i,   0,k) = 0.0d0

        wk(i,   0,k) = coef2*wk(i,   1,k)
        wk(i,  -1,k) = 2.0d0*wk(i,   0,k) - wk(i,   1,k)
        wk(i,  -2,k) = 2.0d0*wk(i,  -1,k) - wk(i,   0,k)
      end do
    end do
!$OMP END PARALLEL DO
  end if

  !================ 上壁 j=nj+1 側 =================
  if (nID(Y_PLUS) .lt. 0) then
!$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,k,phi_av,inv_ls1,inv_ls2,ls_mix,ls_eff,coef1,coef2,x_i,d,w_up,w_dn,w) &
!$OMP& SHARED(ni,nj,nk,dx,dy,uk,vk,wk,phi,uwall,l1,l2,epsb,xup_top,xdn_top,hup_top,hdn_top,tiny,l_free,r0,r1,r2)
    do k=-2, nk+3
      do i=-2, ni+3
        phi_av = 0.5d0*(phi(i,nj,k) + phi(i+1,nj,k))
        inv_ls1 =  phi_av        / max(l1(i,nj,k), tiny)
        inv_ls2 = (1.0d0-phi_av) / max(l2(i,nj,k), tiny)
        ls_mix  = 1.0d0 / max(inv_ls1 + inv_ls2, tiny)

        x_i = (dble(i)+0.5d0)*dx
        w_up = 0.0d0; w_dn = 0.0d0
        if (hup_top(k)) then
          d = x_i - xup_top(k)
          w_up = slip_len_profile(d, r0, r1, r2, l_free, ls_mix)
        end if
        if (hdn_top(k)) then
          d = x_i - xdn_top(k)
          w_dn = slip_len_profile(d, r0, r1, r2, l_free, ls_mix)
        end if

        ls_eff = max(w_up, w_dn)


        ! ---- デバッグ（k = nk/2 のみ、i, x_i, x_up, w）----
        !if (k == nk/2 .and. ls_eff > 0.0d0 .and. hdn_top(k)) then
!$OMP CRITICAL (DBG_TOP_MID)
        !  write(*,'(A,I6,1X,A,F12.5,1X,A,F12.5,1X,A,ES16.8)') &
        !  'i=', i, 'x_i=', x_i, 'x_dn=', xdn_top(k), 'ls=', ls_eff
!$OMP END CRITICAL (DBG_TOP_MID)
        !end if
        !-----------------------------------------------


        coef1 = (2.0d0*dy)         / (2.0d0*ls_eff + dy)
        coef2 = (2.0d0*ls_eff - dy) / (2.0d0*ls_eff + dy)

        uk(i,nj+1,k) = coef1*(+uwall) + coef2*uk(i,nj  ,k)
        uk(i,nj+2,k) = 2.0d0*uk(i,nj+1,k) - uk(i,nj  ,k)
        uk(i,nj+3,k) = 2.0d0*uk(i,nj+2,k) - uk(i,nj+1,k)

        vk(i,nj  ,k) = 0.0d0
        vk(i,nj+1,k) = -vk(i,nj-1,k)
        vk(i,nj+2,k) = -vk(i,nj-2,k)
        vk(i,nj+3,k) = -vk(i,nj-3,k)

        wk(i,nj+1,k) = coef2*wk(i,nj  ,k)
        wk(i,nj+2,k) = 2.0d0*wk(i,nj+1,k) - wk(i,nj  ,k)
        wk(i,nj+3,k) = 2.0d0*wk(i,nj+2,k) - wk(i,nj+1,k)
      end do
    end do
!$OMP END PARALLEL DO
  end if

  !================ x方向の周期コピー =================
!$OMP PARALLEL DO SCHEDULE(static,1) DEFAULT(none) PRIVATE(j,k) SHARED(ni,nj,nk,uk,vk,wk)
  do j = -2, nj+3
    do k = -2, nk+3
      uk(0 ,j,k) = uk(ni+1,j,k)
      uk(-1,j,k) = uk(ni  ,j,k)
      uk(-2,j,k) = uk(ni-1,j,k)
      vk(0 ,j,k) = vk(ni  ,j,k)
      vk(-1,j,k) = vk(ni-1,j,k)
      vk(-2,j,k) = vk(ni-2,j,k)
      wk(0 ,j,k) = wk(ni  ,j,k)
      wk(-1,j,k) = wk(ni-1,j,k)
      wk(-2,j,k) = wk(ni-2,j,k)

      uk(ni+2, j, k) = uk(2, j, k)
      uk(ni+3, j, k) = uk(3, j, k)
      vk(ni+1, j, k) = vk(1, j, k)
      vk(ni+2, j, k) = vk(2, j, k)
      vk(ni+3, j, k) = vk(3, j, k)
      wk(ni+1, j, k) = wk(1, j, k)
      wk(ni+2, j, k) = wk(2, j, k)
      wk(ni+3, j, k) = wk(3, j, k)
    end do
  end do
!$OMP END PARALLEL DO

  !================ z方向の周期コピー =================
!$OMP PARALLEL DO SCHEDULE(static,1) DEFAULT(none) PRIVATE(i,j) SHARED(ni,nj,nk,uk,vk,wk)
  do i = -2, ni+3
    do j = -2, nj+3
      uk(i, j,   0) = uk(i, j, nk  )
      uk(i, j,  -1) = uk(i, j, nk-1)
      uk(i, j,  -2) = uk(i, j, nk-2)
      vk(i, j,   0) = vk(i, j, nk  )
      vk(i, j,  -1) = vk(i, j, nk-1)
      vk(i, j,  -2) = vk(i, j, nk-2)
      wk(i, j,   0) = wk(i, j, nk+1)
      wk(i, j,  -1) = wk(i, j, nk  )
      wk(i, j,  -2) = wk(i, j, nk-1)

      uk(i, j, nk+1) = uk(i, j,   1)
      uk(i, j, nk+2) = uk(i, j,   2)
      uk(i, j, nk+3) = uk(i, j,   3)
      vk(i, j, nk+1) = vk(i, j,   1)
      vk(i, j, nk+2) = vk(i, j,   2)
      vk(i, j, nk+3) = vk(i, j,   3)
      wk(i, j, nk+2) = wk(i, j,   2)
      wk(i, j, nk+3) = wk(i, j,   3)
    end do
  end do
!$OMP END PARALLEL DO

end subroutine bndu_localized_free