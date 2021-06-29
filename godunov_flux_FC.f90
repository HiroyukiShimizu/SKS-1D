subroutine godunov_flux_FC(h_L,rho_L,u_L,na_L,ns_L,nw_L,aaa_L,Fr_L,enthal_L,&
                           f1,f2,f3,f4,f5,f6,theta)

  use parameter
  implicit none

  double precision, intent(in) :: h_L,rho_L,u_L,na_L,ns_L,nw_L,aaa_L,Fr_L,enthal_L
  double precision, intent(in) :: theta
  double precision, intent(out) :: f1,f2,f3,f4,f5,f6
  double precision :: h,rho,u,na,ns,nw,enthal
  double precision :: h_nondim, u_nondim
  double precision :: h_N_nondim, u_N_nondim, Fr_N

  rho = rho_L
  na = na_L
  ns = ns_L
  nw = nw_L
  enthal = enthal_L

  Fr_N = Fr_N0 * dsqrt(rho/rho_a)

  if(Fr_L<=Fr_N) then
    if(1.d0<=Fr_L .and. Fr_L<=Fr_N) then
      h = h_L
      u = u_L
    elseif(Fr_L<1.d0 .and. 1.d0<Fr_N) then
      h_nondim = (Fr_L+2.d0)**2 / 9.d0
      u_nondim = (Fr_L+2.d0) / 3.d0
      h = h_L * h_nondim
      u = aaa_L * u_nondim
    else
      h_N_nondim = ( (2.d0+Fr_L)/(2.d0+Fr_N) )**2
      u_N_nondim = Fr_N * (2.d0+Fr_L)/(2.d0+Fr_N)
      h = h_L * h_N_nondim
      u = aaa_L * u_N_nondim
    endif
  else
    h = h_L
    u = u_L
  endif

  f1 = rho * u * h
  f2 = na * rho * u * h
  f3 = ns * rho * u * h
  f4 = rho*u*u*h + 0.5d0*(rho-rho_a)*grav*h*h*dcos(theta)
  f5 = rho * u * enthal * h
  f6 = nw * rho * u * h


  return
end subroutine godunov_flux_FC
