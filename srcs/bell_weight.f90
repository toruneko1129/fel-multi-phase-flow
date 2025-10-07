pure elemental real*8 function bell_weight(d, epsb)
  implicit none
  real*8, intent(in) :: d, epsb
  real*8 :: r, pi
  pi = acos(-1.0d0)

  if (epsb <= 0.0d0) then
    bell_weight = 0.0d0
  else if (abs(d) >= epsb) then
    bell_weight = 0.0d0
  else
    r = d/epsb
    bell_weight = 0.5d0*(1.0d0 + cos(pi*r))
    bell_weight = bell_weight*bell_weight
  end if
end function bell_weight