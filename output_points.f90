subroutine output_points(n,time,i_output_1,i_output_2,i_output_3,i_output_4)

  use parameter
  use variable
  implicit none

  integer, intent(in) :: n,i_output_1,i_output_2,i_output_3,i_output_4
  double precision, intent(in) :: time
  double precision :: h_1, rho_1, u_1, n_a1, n_s1, phi_a1, phi_s1, T_1, Fr_1, Ri_1
  double precision :: hH_1, uH_1, FrH_1, RiH_1, theta_1, zc_1, zb_1, Sa_1
  double precision :: h_2, rho_2, u_2, n_a2, n_s2, phi_a2, phi_s2, T_2, Fr_2, Ri_2
  double precision :: hH_2, uH_2, FrH_2, RiH_2, theta_2, zc_2, zb_2, Sa_2
  double precision :: h_3, rho_3, u_3, n_a3, n_s3, phi_a3, phi_s3, T_3, Fr_3, Ri_3
  double precision :: hH_3, uH_3, FrH_3, RiH_3, theta_3, zc_3, zb_3, Sa_3
  double precision :: h_4, rho_4, u_4, n_a4, n_s4, phi_a4, phi_s4, T_4, Fr_4, Ri_4
  double precision :: hH_4, uH_4, FrH_4, RiH_4, theta_4, zc_4, zb_4, Sa_4


  if(lift_off==0) then
    !*** dilute (L) ***
    if(time==0.d0) then
      write(600,*)'# x=',x(i_output_1),'    x-x_0=',x(i_output_1)-x_0
      write(600,*)'# t, h, rho, u, n_a, n_s, T, Fr, Ri, rho/rho_a, phi_a, phi_s'
      write(700,*)'# x=',x(i_output_2),'    x-x_0=',x(i_output_2)-x_0
      write(700,*)'# t, h, rho, u, n_a, n_s, T, Fr, Ri, rho/rho_a, phi_a, phi_s'
      write(800,*)'# x=',x(i_output_3),'    x-x_0=',x(i_output_3)-x_0
      write(800,*)'# t, h, rho, u, n_a, n_s, T, Fr, Ri, rho/rho_a, phi_a, phi_s'
      write(900,*)'# x=',x(i_output_4),'    x-x_0=',x(i_output_4)-x_0
      write(900,*)'# t, h, rho, u, n_a, n_s, T, Fr, Ri, rho/rho_a, phi_a, phi_s'
    endif
    h_1 = h(i_output_1)
    rho_1 = rho(i_output_1)
    u_1 = u(i_output_1)
    n_a1 = na(i_output_1)
    n_s1 = ns(i_output_1)
    phi_a1 = phia(i_output_1)
    phi_s1 = phis(i_output_1)
    T_1 = T(i_output_1)
    Fr_1 = Fr(i_output_1)
    if(Fr_1>0.d0) then
      Ri_1 = 1.d0 / Fr_1 / Fr_1
    else
      Ri_1 = 0.d0
    endif

    h_2 = h(i_output_2)
    rho_2 = rho(i_output_2)
    u_2 = u(i_output_2)
    n_a2 = na(i_output_2)
    n_s2 = ns(i_output_2)
    phi_a2 = phia(i_output_2)
    phi_s2 = phis(i_output_2)
    T_2 = T(i_output_2)
    Fr_2 = Fr(i_output_2)
    if(Fr_2>0.d0) then
      Ri_2 = 1.d0 / Fr_2 / Fr_2
    else
      Ri_2 = 0.d0
    endif

    h_3 = h(i_output_3)
    rho_3 = rho(i_output_3)
    u_3 = u(i_output_3)
    n_a3 = na(i_output_3)
    n_s3 = ns(i_output_3)
    phi_a3 = phia(i_output_3)
    phi_s3 = phis(i_output_3)
    T_3 = T(i_output_3)
    Fr_3 = Fr(i_output_3)
    if(Fr_3>0.d0) then
      Ri_3 = 1.d0 / Fr_3 / Fr_3
    else
      Ri_3 = 0.d0
    endif

    h_4 = h(i_output_4)
    rho_4 = rho(i_output_4)
    u_4 = u(i_output_4)
    n_a4 = na(i_output_4)
    n_s4 = ns(i_output_4)
    phi_a4 = phia(i_output_4)
    phi_s4 = phis(i_output_4)
    T_4 = T(i_output_4)
    Fr_4 = Fr(i_output_4)
    if(Fr_4>0.d0) then
      Ri_4 = 1.d0 / Fr_4 / Fr_4
    else
      Ri_4 = 0.d0
    endif

    if(non_dim_output==0) then
      write(600,'(12E26.16)') time, h_1, rho_1, u_1, n_a1, n_s1, T_1, &
                             Fr_1, Ri_1, rho_1/rho_a, phi_a1, phi_s1
      write(700,'(12E26.16)') time, h_2, rho_2, u_2, n_a2, n_s2, T_2, &
                             Fr_2, Ri_2, rho_2/rho_a, phi_a2, phi_s2
      write(800,'(12E26.16)') time, h_3, rho_3, u_3, n_a3, n_s3, T_3, &
                             Fr_3, Ri_3, rho_3/rho_a, phi_a3, phi_s3
      write(900,'(12E26.16)') time, h_4, rho_4, u_4, n_a4, n_s4, T_4, &
                             Fr_4, Ri_4, rho_4/rho_a, phi_a4, phi_s4

    elseif(non_dim_output==1) then
      write(600,'(12E26.16)') time/T_char, h_1/H_char, rho_1/rho_0, &
                             u_1/U_char, n_a1, n_s1, T_1/T_0, &
                             Fr_1, Ri_1, rho_1/rho_a, phi_a1, phi_s1
      write(700,'(12E26.16)') time/T_char, h_2/H_char, rho_2/rho_0, &
                             u_2/U_char, n_a2, n_s2, T_2/T_0, &
                             Fr_2, Ri_2, rho_2/rho_a, phi_a2, phi_s2
      write(800,'(12E26.16)') time/T_char, h_3/H_char, rho_3/rho_0, &
                             u_3/U_char, n_a3, n_s3, T_3/T_0, &
                             Fr_3, Ri_3, rho_3/rho_a, phi_a3, phi_s3
      write(900,'(12E26.16)') time/T_char, h_4/H_char, rho_4/rho_0, &
                             u_4/U_char, n_a4, n_s4, T_4/T_0, &
                             Fr_4, Ri_4, rho_4/rho_a, phi_a4, phi_s4

    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif

  endif

  !*** dense (H) & deposit (D) ***
  if(n==0) then
    write(650,*)'# x=',x(i_output_1),'    x-x_0=',x(i_output_1)-x_0
    write(650,*)'# t, h_H, u_H, Fr_H, Ri_H, z_c, Sa, z_b'
    write(750,*)'# x=',x(i_output_2),'    x-x_0=',x(i_output_2)-x_0
    write(750,*)'# t, h_H, u_H, Fr_H, Ri_H, z_c, Sa, z_b'
    write(850,*)'# x=',x(i_output_3),'    x-x_0=',x(i_output_3)-x_0
    write(850,*)'# t, h_H, u_H, Fr_H, Ri_H, z_c, Sa, z_b'
    write(950,*)'# x=',x(i_output_4),'    x-x_0=',x(i_output_4)-x_0
    write(950,*)'# t, h_H, u_H, Fr_H, Ri_H, z_c, Sa, z_b'
  endif

  hH_1 = hH(i_output_1)
  uH_1 = uH(i_output_1)
  FrH_1 = FrH(i_output_1)
  if(FrH_1<=0.d0) RiH_1 = 0.d0
  if(FrH_1>0.d0) RiH_1 = 1.d0 / FrH_1 / FrH_1
  theta_1 = theta(i_output_1)
  if(Ws_type==1.or.Ws_type==2) then
    Sa_1 = rho_s * (diameter*uH_1/hH_1)**2 / ( (rho_s-rho_a)*grav*hH_1*dcos(theta_1) )
  else
    Sa_1 = 0.d0
  endif
  zc_1 = z_c(i_output_1)
  zb_1 = z_b(i_output_1)

  hH_2 = hH(i_output_2)
  uH_2 = uH(i_output_2)
  FrH_2 = FrH(i_output_2)
  if(FrH_2<=0.d0) RiH_2 = 0.d0
  if(FrH_2>0.d0) RiH_2 = 1.d0 / FrH_2 / FrH_2
  theta_2 = theta(i_output_2)
  if(Ws_type==1.or.Ws_type==2) then
    Sa_2 = rho_s * (diameter*uH_2/hH_2)**2 / ( (rho_s-rho_a)*grav*hH_2*dcos(theta_2) )
  else
    Sa_2 = 0.d0
  endif
  zc_2 = z_c(i_output_2)
  zb_2 = z_b(i_output_2)

  hH_3 = hH(i_output_3)
  uH_3 = uH(i_output_3)
  FrH_3 = FrH(i_output_3)
  if(FrH_3<=0.d0) RiH_3 = 0.d0
  if(FrH_3>0.d0) RiH_3 = 1.d0 / FrH_3 / FrH_3
  theta_3 = theta(i_output_3)
  if(Ws_type==1.or.Ws_type==2) then
    Sa_3 = rho_s * (diameter*uH_3/hH_3)**2 / ( (rho_s-rho_a)*grav*hH_3*dcos(theta_3) )
  else
    Sa_3 = 0.d0
  endif
  zc_3 = z_c(i_output_3)
  zb_3 = z_b(i_output_3)

  hH_4 = hH(i_output_4)
  uH_4 = uH(i_output_4)
  FrH_4 = FrH(i_output_4)
  if(FrH_4<=0.d0) RiH_4 = 0.d0
  if(FrH_4>0.d0) RiH_4 = 1.d0 / FrH_4 / FrH_4
  theta_4 = theta(i_output_4)
  if(Ws_type==1.or.Ws_type==2) then
    Sa_4 = rho_s * (diameter*uH_4/hH_4)**2 / ( (rho_s-rho_a)*grav*hH_4*dcos(theta_4) )
  else
    Sa_4 = 0.d0
  endif
  zc_4 = z_c(i_output_4)
  zb_4 = z_b(i_output_4)

  if(non_dim_output==0) then
    write(650,'(8E26.16)') time, hH_1, uH_1, FrH_1, RiH_1, zc_1, Sa_1, zb_1
    write(750,'(8E26.16)') time, hH_2, uH_2, FrH_2, RiH_2, zc_2, Sa_2, zb_2
    write(850,'(8E26.16)') time, hH_3, uH_3, FrH_3, RiH_3, zc_3, Sa_3, zb_3
    write(950,'(8E26.16)') time, hH_4, uH_4, FrH_4, RiH_4, zc_4, Sa_4, zb_4
  elseif(non_dim_output==1) then
    write(650,'(8E26.16)') time/T_char, hH_1/H_char, uH_1/U_char, FrH_1, &
                           RiH_1, zc_1/H_char, Sa_1, zb_1/H_char
    write(750,'(8E26.16)') time/T_char, hH_2/H_char, uH_2/U_char, FrH_2, &
                           RiH_2, zc_2/H_char, Sa_2, zb_2/H_char
    write(850,'(8E26.16)') time/T_char, hH_3/H_char, uH_3/U_char, FrH_3, &
                           RiH_3, zc_3/H_char, Sa_3, zb_3/H_char
    write(950,'(8E26.16)') time/T_char, hH_4/H_char, uH_4/U_char, FrH_4, &
                           RiH_4, zc_4/H_char, Sa_4, zb_4/H_char
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif


  return
end subroutine output_points
