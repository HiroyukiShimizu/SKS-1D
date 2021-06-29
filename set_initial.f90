subroutine set_initial()

  use parameter
  use variable
  implicit none

  integer :: i

  !*** dilute (L) ***
  if(rho_0/rho_a<1.d0) then
    write(*,*)'*** ERROR **** (LIFT OFF) in set_initial.f90'
    write(*,*)'rho_0/rho_a=',rho_0/rho_a
    stop
  else
    lift_off = 0.d0
  endif

  lift_off_FC = 0
  FC = msource + 2
  x_N_steady = 0.d0 + left_boundary_point
  i_front_steady = 0
  dx_FC = 0.d0
  x_N = dble(FC-1)*dx + dx_FC + left_boundary_point
  x_N_max = x_N
  i_front_L = FC-1

  do i = i_start-1, i_front_L
    h(i) = h_0
    rho(i) = rho_0
    u(i) = u_0
    na(i) = n_a0
    ns(i) = n_s0
    nw(i) = n_w0
    phia(i) = phi_a0
    phis(i) = phi_s0
    phiw(i) = phi_w0
    aaa(i) = dsqrt(h_0*grav*dcos(theta(i))*(rho_0-rho_a)/rho_0)
    Fr(i) = dabs(u_0) / aaa(i)
    Cp(i) = C_p0
    Cp_old(i) = C_p0
    T(i) = T_0
    T_boil(i) = T_boil0
    c(i) = 0.d0
    if(thermal_type==0) then
      enthal(i) = 0.d0
    elseif(thermal_type/=0) then
      enthal(i) = C_p0 * T_0
    endif

    Q(i,1) = rho(i) * h(i)
    Q(i,2) = na(i) * rho(i) * h(i)
    Q(i,3) = ns(i) * rho(i) * h(i)
    Q(i,4) = rho(i) * u(i) * h(i)
    Q(i,6) = nw(i) * rho(i) * h(i)
    if(thermal_type==0) then
      Q(i,5) = 0.d0
    elseif(thermal_type/=0) then
      Q(i,5) = rho(i) * Cp(i) * T(i) * h(i)
    endif
  enddo

  do i = i_front_L+1, mx+1
    h(i) = 0.d0
    rho(i) = 0.d0
    u(i) = 0.d0
    na(i) = 0.d0
    ns(i) = 0.d0
    nw(i) = 0.d0
    phia(i) = 0.d0
    phis(i) = 0.d0
    phiw(i) = 0.d0
    aaa(i) = 0.d0
    Fr(i) = 0.d0
    T(i) = 0.d0
    T_boil(i) = 0.d0
    enthal(i) = 0.d0
    Cp(i) = 0.d0
    Cp_old(i) = 0.d0
    c(i) = 0.d0
    Q(i,1) = 0.d0
    Q(i,2) = 0.d0
    Q(i,3) = 0.d0
    Q(i,4) = 0.d0
    Q(i,5) = 0.d0
    Q(i,6) = 0.d0
  enddo

  !*** dense ***
  if(h_H0>0.d0) then
    i_front_H = i_front_L
    do i = i_start-1, i_front_H
      hH(i) = h_H0
      uH(i) = u_H0
      aaaH(i) = dsqrt(h_H0*g_H*dcos(theta(i)))
      FrH(i) = u_H0/aaaH(i)
      QH(i,1) = h_H0 + h_Hchar*eps ! + artificial bed
      QH(i,2) = (h_H0+h_Hchar*eps)*u_H0 ! + artificial bed
    enddo
    do i = i_front_H+1, mx+1
      hH(i) = 0.d0
      uH(i) = 0.d0
      aaaH(i) = 0.d0
      FrH(i) = 0.d0
      QH(i,1) = eps*h_Hchar ! artificial bed
      QH(i,2) = 0.d0
    enddo

  else
    i_front_H = i_start - 1
    hH(:) = 0.d0
    uH(:) = 0.d0
    aaaH(:) = 0.d0
    FrH(:) = 0.d0
    QH(:,1) = eps*h_Hchar ! artificial bed
    QH(:,2) = 0.d0
  endif
  x_NH = x(i_front_H)
  x_NH_max = x_NH
  x_NH_steady = left_boundary_point
  i_front_H_steady = 0

  !*** Deposit and topography ***
  i_front_D = i_start - 1
  x_ND = x(i_front_D)
  hD(:) = 0.d0
  z_b0(:) = 0.d0

  z_b(:) = z_b0(:) + hD(:)
  z_c(:) = hH(:) + z_b(:)


  return
end subroutine set_initial
