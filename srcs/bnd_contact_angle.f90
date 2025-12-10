subroutine bnd_contact_angle(nID, ni, nj, nk, qk, theta_array, dx, dy, dz)
  implicit none
  include 'param.h'

  integer, intent(in)    :: nID(6), ni, nj, nk
  real*8, intent(inout)  :: qk(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)     :: theta_array(-2:ni+3, -2:nj+3, -2:nk+3)
  real*8, intent(in)     :: dx, dy, dz

  integer :: i, k
  real*8 :: pi, eps, theta_rad, s, c, cot_t
  real*8 :: gx, gz, gt, dqdY

  pi  = acos(-1.0d0)
  eps = 1.0d-14   ! 1e-20 でもいいが、勾配計算だとこれくらいの方が安定しやすい

  !======================================================
  ! 下壁 (Y_MINUS): 壁法線 = +y (s_w = +1)
  ! 勾配は interior の j=1 で評価
  !======================================================
  if (nID(Y_MINUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,s,c,cot_t,gx,gz,gt,dqdY) &
!$OMP& SHARED(ni,nk,qk,theta_array,dx,dy,dz,pi,eps)
    do k = 1, nk
      do i = 1, ni
        theta_rad = theta_array(i,1,k) * (pi/180.0d0)
        s = sin(theta_rad)
        c = cos(theta_rad)

        ! 壁面内（x,z）の接線勾配：中央差分（周期ゴーストが先に入っている前提）
        gx = (qk(i+1,1,k)   - qk(i-1,1,k))   / (2.0d0*dx)
        !gz = (qk(i,  1,k+1) - qk(i,  1,k-1)) / (2.0d0*dz)
        gz = 0.0d0
        gt = sqrt(gx*gx + gz*gz)

        ! sin(theta)->0 や gt->0 のときは爆発するので、自然に Neumann(∂y=0)へ退避
        if (gt .le. eps .or. abs(s) .le. eps) then
          qk(i,0,k)  = qk(i,1,k)
        else
          cot_t = c / s
          ! dqdY = -s_w*gt*cot(theta),  下壁は s_w=+1
          dqdY = - gt * cot_t
          ! (q1 - q0)/dy = dqdY  ->  q0 = q1 - dy*dqdY
          qk(i,0,k)  = qk(i,1,k) - dy*dqdY
        end if

        ! 追加ゴースト：線形外挿（必要な段数だけ）
        qk(i,-1,k) = 2.0d0*qk(i,0,k)  - qk(i,1,k)
        qk(i,-2,k) = 2.0d0*qk(i,-1,k) - qk(i,0,k)

        ! VOFなら念のためクリップ（任意）
        ! qk(i,0,k)  = max(0.0d0, min(1.0d0, qk(i,0,k)))
        ! qk(i,-1,k) = max(0.0d0, min(1.0d0, qk(i,-1,k)))
        ! qk(i,-2,k) = max(0.0d0, min(1.0d0, qk(i,-2,k)))
      end do
    end do
!$OMP END PARALLEL DO
  end if


  !======================================================
  ! 上壁 (Y_PLUS): 壁法線 = -y (s_w = -1)
  ! 勾配は interior の j=nj で評価
  ! ★あなたの現行コードのバグ：上壁なのに nq_y=cosθ のまま（符号が違う）
  !======================================================
  if (nID(Y_PLUS) .lt. 0) then
!$OMP PARALLEL DO PRIVATE(i,k,theta_rad,s,c,cot_t,gx,gz,gt,dqdY) &
!$OMP& SHARED(ni,nk,nj,qk,theta_array,dx,dy,dz,pi,eps)
    do k = 1, nk
      do i = 1, ni
        theta_rad = theta_array(i,nj,k) * (pi/180.0d0)
        s = sin(theta_rad)
        c = cos(theta_rad)

        gx = (qk(i+1,nj,k)   - qk(i-1,nj,k))   / (2.0d0*dx)
        !gz = (qk(i,  nj,k+1) - qk(i,  nj,k-1)) / (2.0d0*dz)
        gz = 0.0d0
        gt = sqrt(gx*gx + gz*gz)

        if (gt .le. eps .or. abs(s) .le. eps) then
          qk(i,nj+1,k) = qk(i,nj,k)
        else
          cot_t = c / s
          ! dqdY = -s_w*gt*cot(theta),  上壁は s_w=-1 -> dqdY = +gt*cot
          dqdY = + gt * cot_t
          ! (q_{nj+1}-q_{nj})/dy = dqdY -> q_{nj+1} = q_{nj} + dy*dqdY
          qk(i,nj+1,k) = qk(i,nj,k) + dy*dqdY
        end if

        qk(i,nj+2,k) = 2.0d0*qk(i,nj+1,k) - qk(i,nj,k)
        qk(i,nj+3,k) = 2.0d0*qk(i,nj+2,k) - qk(i,nj+1,k)

        ! VOFなら任意でクリップ
        ! qk(i,nj+1,k) = max(0.0d0, min(1.0d0, qk(i,nj+1,k)))
        ! qk(i,nj+2,k) = max(0.0d0, min(1.0d0, qk(i,nj+2,k)))
        ! qk(i,nj+3,k) = max(0.0d0, min(1.0d0, qk(i,nj+3,k)))
      end do
    end do
!$OMP END PARALLEL DO
  end if

  return
end subroutine bnd_contact_angle
