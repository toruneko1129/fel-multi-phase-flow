pure elemental real*8 function slip_len_profile(d, r0, r1, r2, l_free, l_mix) result(ls_eff)
  implicit none
  real*8, intent(in) :: d, r0, r1, r2, l_free, l_mix
  real*8 :: a, t, pi
  real*8 :: eps

  pi = acos(-1.0d0)
  eps = 1.0d-12

  ! パラメータチェック
  if (r0 <= 0.0d0 .or. r1 <= r0 .or. r2 <= r1) then
    ls_eff = l_mix
    return
  end if

  a = abs(d)

  if (a <= r0 + eps) then
    ! コア：free-slip 相当
    ls_eff = l_free

  elseif (a < r1 + eps) then
    ! コア → mix への遷移（C^1 連続の余弦窓）
    t = (a - r0) / (r1 - r0)
    t = 0.5d0*(1.0d0 + cos(pi*t))
    t = t*t
    ls_eff = t*l_free + (1.0d0 - t)*l_mix

  elseif (a < r2 + eps) then
    ! mix → no-slip への遷移
    t = (a - r1) / (r2 - r1)
    t = 0.5d0*(1.0d0 + cos(pi*t))
    t = t*t
    ls_eff = t*l_mix                    ! + (1-t)*0

  else
    ! 遠方：no-slip
    ls_eff = 0.0d0
  end if
end function slip_len_profile
