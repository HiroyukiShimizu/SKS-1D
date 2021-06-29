subroutine cons_to_prim_FC(A1, A2, A3, A4, A5, A6, dt, &
                           h_FC, rho_FC, u_FC, &
                           na_FC, ns_FC, nw_FC, phia_FC, phis_FC, phiw_FC, &
                           T_FC, enthal_FC, Cp_FC, aaa_FC, Fr_FC, &
                           theta_FC, Tboil_FC)

  use parameter
  use variable
  implicit none

  double precision, intent(in) :: A1,A2,A3,A4,A5,A6,dt,theta_FC
  double precision, intent(out) :: h_FC,rho_FC,u_FC
  double precision, intent(out) :: na_FC,ns_FC,nw_FC,phia_FC,phis_FC,phiw_FC
  double precision, intent(out) :: T_FC,enthal_FC,Cp_FC,aaa_FC,Fr_FC
  double precision, intent(out) :: Tboil_FC

  double precision :: Fr_N
  double precision :: pv_FC
  double precision, dimension(1:3) :: WW ! = (h_FC,u_FC,dx_FC_new)

  na_FC = A2 / A1
  ns_FC = A3 / A1
  nw_FC = A6 / A1
  Cp_FC = ns_FC*C_s + nw_FC*C_pw + na_FC*C_pa + (1.d0-ns_FC-nw_FC-na_FC)*C_pv
  if(thermal_type==0) then
    enthal_FC = 0.d0
    T_FC = T_0
  elseif(thermal_type==1) then
    enthal_FC = A5 / A1
    T_FC = T_star(FC) - dt/dx*u_star(FC)*(T_star(FC)-T_star(FC-1)) &
           +(enthal_FC-enthal_star(FC))/Cp_FC &
           +dt/dx*u_star(FC)*(enthal_star(FC)-enthal_star(FC-1))/Cp_FC
  elseif(thermal_type==2) then
    enthal_FC = A5 / A1
    T_FC = enthal_FC/Cp_FC
  endif


  pv_FC = p_a0 * (1.d0-ns_FC-nw_FC-na_FC)*R_v
  pv_FC = pv_FC / ( (1.d0-ns_FC-nw_FC-na_FC)*R_v + na_FC*R_a )
  if(1.d0-ns_FC-nw_FC-na_FC <= 0.d0) then
    pv_FC = 0.d0
    Tboil_FC = 0.d0
  endif


  rho_FC = ns_FC/rho_s + nw_FC/rho_w + ( na_FC*R_a + (1.d0-ns_FC-nw_FC-na_FC)*R_v ) * T_FC / p_a0
  rho_FC = 1.d0 / rho_FC

  if(dx_FC<=0.d0.and.rho_FC<=rho_a*(1.d0+1.d-2)) then
    lift_off_FC = 1
    if(R_v==R_a.and.thermal_type==0) then
      write(*,*)'rho_FC=',rho_FC
      stop
    endif
    return
  endif
  if(dx_FC>0.d0.and.rho_FC<=rho_a*(1.d0+1.d-3)) then
    lift_off_FC = 1
    if(R_v==R_a.and.thermal_type==0) then
      write(*,*)'rho_FC=',rho_FC
      stop
    endif
    return
  endif


  lift_off_FC = 0

  phia_FC = na_FC*rho_FC * R_a*T_FC/p_a0
  phis_FC = ns_FC*rho_FC / rho_s
  phiw_FC = nw_FC*rho_FC / rho_w

  !*** initial condition of Newton-Raphson method ***
  Fr_N = Fr_N0*dsqrt(rho_FC/rho_a)
  if(FC_type==3) Fr_N=1.d0
  if(dx_FC<=0.d0) then
    if( Fr_star(FC-1) <= Fr_N ) then
      WW(1) = h_star(FC-1) * ( (2.d0+Fr_star(FC-1)) / (2.d0+Fr_N) )**2
    else
      WW(1) = h_star(FC-1) !!!! tentative
    endif
    WW(2) = Fr_N * dsqrt( WW(1) * grav * dcos(theta_FC) * (rho_FC-rho_a)/rho_FC )
    WW(3) = A1 / (rho_FC*WW(1))
  else
    WW(1) = h_star(FC)
    WW(2) = u_star(FC)
    WW(3) = A1 / (rho_FC*WW(1))
  endif

    if(FC_type==1) then
      call not_kinematic_condition(WW,rho_FC,dt,A1,A4,theta_FC)
    elseif(FC_type==2) then
      call not_momentum_equation(WW,rho_FC,dt,A1,dx_FC,theta_FC)
    elseif(FC_type==3) then
      call not_kinematic_condition(WW,rho_FC,dt,A1,A4,theta_FC)
      if(WW(3)<=dx_FC) then
        WW(1) = h_star(FC)
        WW(2) = u_star(FC)
        WW(3) = A1 / (rho_FC*WW(1))
        call not_momentum_equation(WW,rho_FC,dt,A1,dx_FC,theta_FC)
      endif
    elseif(FC_type==4) then
      call not_front_condition(WW,rho_FC,dt,A1,A4,dx_FC,theta_FC)
    else
      write(*,*)'*** ERROR *** (FC_type)'
      write(*,*)'  FC_type=',FC_type
      stop
    endif

  h_FC = WW(1)
  u_FC = WW(2)
  dx_FC_new = WW(3)

  if(dx_FC_new<=dx_FC) then
    dx_FC_new = 0.d0
    h_FC = 0.d0
    u_FC = 0.d0
    lift_off_FC = 1
  endif

  aaa_FC = dsqrt(h_FC*grav*dcos(theta_FC)*(rho_FC-rho_a)/rho_FC)
  Fr_FC = dabs(u_FC) / aaa_FC


  return
end subroutine cons_to_prim_FC
