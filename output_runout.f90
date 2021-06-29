subroutine output_runout()

  use parameter
  use variable
  implicit none

  integer :: i
  double precision :: u_ave, E_ave, Eu_ave

  u_ave = 0.d0
  E_ave = 0.d0
  Eu_ave = 0.d0

  do i = i_start, i_front_L
    u_ave = u_ave + dabs(u(i))
    E_ave = E_ave + Ent_func(Fr(i))
    if(configuration_type == 1) then
      Eu_ave = Eu_ave + Ent_func(Fr(i))*dabs(u(i))
    elseif(configuration_type == 2) then
      Eu_ave = Eu_ave + Ent_func(Fr(i))*dabs(u(i))*x(i)*dx
    endif
  enddo
  u_ave = u_ave / dble(i_front_L-i_start+1)
  E_ave = E_ave / dble(i_front_L-i_start+1)
  if(configuration_type == 1) then
    Eu_ave = Eu_ave / dble(i_front_L-i_start+1)
  elseif(configuration_type == 2) then
    Eu_ave = Eu_ave * 2.d0 / (x_N_steady*x_N_steady-x_0*x_0)
  endif

  open(unit=100,file='runout.nml',status='unknown',form='formatted')

  if(non_dim_output==0) then
    write(100,*)'&dilute_run_list'
    write(100,*)' x_N_max = ', x_N_max,','
    write(100,*)' x_N_steady = ', x_N_steady,','
    write(100,*)' h_N_steady = ', h_N_steady,','
    write(100,*)' rho_N_steady = ', rho_N_steady,','
    write(100,*)' u_N_steady = ', u_N_steady,','
    write(100,*)' na_N_steady = ', na_N_steady,','
    write(100,*)' ns_N_steady = ', ns_N_steady,','
    write(100,*)' T_N_steady = ', T_N_steady,','
    write(100,*)' Fr_N_steady = ', Fr_N_steady,','
    write(100,*)' Ri_N_steady = ', 1.d0/(Fr_N_steady*Fr_N_steady),','
    write(100,*)' enthal_N_steady = ', enthal_N_steady,','
    write(100,*)' phia_N_steady = ', phia_N_steady,','
    write(100,*)' phis_N_steady = ', phis_N_steady,','
    write(100,*)' time_steady = ', time_steady,','
    write(100,*)' T_N_steady_T_a = ', T_N_steady/T_a,','
    write(100,*)' M_dot_top = ', M_dot_top,','
    write(100,*)' M_dot_base = ', M_dot_base
    write(100,*)'/'
    write(100,*)''
    write(100,*)'&dense_run_list'
    write(100,*)' x_NH_max = ', x_NH_max,','
    write(100,*)' x_NH_steady = ', x_NH_steady,','
    write(100,*)' time_H_steady = ', time_H_steady,','
    write(100,*)' M_dot_top_H = ', M_dot_top_H,','
    write(100,*)' M_dot_base_H = ', M_dot_base_H
    write(100,*)'/'
    write(100,*)''
    write(100,*)'&dilute_ent_list'
    write(100,*)' u_ave = ', u_ave,','
    write(100,*)' E_ave = ', E_ave,','
    write(100,*)' Eu_ave = ', Eu_ave
    write(100,*)'/'
    write(100,*)''

  elseif(non_dim_output==1) then
    write(100,*)'&dilute_run_list'
    write(100,*)' x_N_max = ', x_N_max / x_0,','
    write(100,*)' x_N_steady = ', x_N_steady / x_0,','
    write(100,*)' h_N_steady = ', h_N_steady / H_char,','
    write(100,*)' rho_N_steady = ', rho_N_steady / rho_0,','
    write(100,*)' u_N_steady = ', u_N_steady / U_char,','
    write(100,*)' na_N_steady = ', na_N_steady,','
    write(100,*)' ns_N_steady = ', ns_N_steady,','
    write(100,*)' T_N_steady = ', T_N_steady / T_0,','
    write(100,*)' Fr_N_steady = ', Fr_N_steady,','
    write(100,*)' Ri_N_steady = ', 1.d0/(Fr_N_steady*Fr_N_steady),','
    write(100,*)' enthal_N_steady = ', enthal_N_steady / (C_p0*T_0),','
    write(100,*)' phia_N_steady = ', phia_N_steady,','
    write(100,*)' phis_N_steady = ', phis_N_steady,','
    write(100,*)' time_steady = ', time_steady / T_char,','
    write(100,*)' T_N_steady_T_a = ', T_N_steady / T_a,','
    write(100,*)' M_dot_top = ', M_dot_top / M_dot_0,','
    write(100,*)' M_dot_base = ', M_dot_base / M_dot_0
    write(100,*)'/'
    write(100,*)''
    write(100,*)'&dense_run_list'
    write(100,*)' x_NH_max = ', x_NH_max / x_0,','
    write(100,*)' x_NH_steady = ', x_NH_steady / x_0,','
    write(100,*)' time_H_steady = ', time_H_steady / T_char,','
    write(100,*)' M_dot_top_H = ', M_dot_top_H / M_dot_0,','
    write(100,*)' M_dot_base_H = ', M_dot_base_H / M_dot_0
    write(100,*)'/'
    write(100,*)''
    write(100,*)''
    write(100,*)'&dilute_ent_list'
    write(100,*)' u_ave = ', u_ave / U_char,','
    write(100,*)' E_ave = ', E_ave,','
    write(100,*)' Eu_ave = ', Eu_ave / U_char
    write(100,*)'/'
    write(100,*)''

  else
    write(*,*)'*** ERROR *** (in output_runout)'
    write(*,*)'   non_dim_output =', non_dim_output
    stop
  endif

  close(100)

  return
end subroutine output_runout
