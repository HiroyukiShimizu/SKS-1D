subroutine output_front(n,time)

  use parameter
  use variable
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: time
  double precision :: h_FC, rho_FC, u_FC, n_aFC, n_sFC, phi_aFC, phi_sFC, T_FC, Fr_FC, Ri_FC


  if(lift_off==0) then
    !*** dilute (L) ***
    if(time==0.d0) then
      write(60,*)'# t, x_N, h_FC, rho_FC, u_FC, n_aFC, n_sFC, T_FC, Fr_FC, Ri_FC, rho_FC/rho_a, phi_aFC, phi_sFC, x_N-x_0'
    endif
    h_FC = h(i_front_L)
    rho_FC = rho(i_front_L)
    u_FC = u(i_front_L)
    n_aFC = na(i_front_L)
    n_sFC = ns(i_front_L)
    phi_aFC = phia(i_front_L)
    phi_sFC = phis(i_front_L)
    T_FC = T(i_front_L)
    Fr_FC = Fr(i_front_L)
    if(Fr_FC>0.d0) then
      Ri_FC = 1.d0 / Fr_FC / Fr_FC
    else
      Ri_FC = 0.d0
    endif
    if(non_dim_output==0) then
      write(60,'(14E26.16)') time, x_N, h_FC, rho_FC, u_FC, n_aFC, n_sFC, T_FC, &
                             Fr_FC, Ri_FC, rho_FC/rho_a, phi_aFC, phi_sFC, x_N-x_0
    elseif(non_dim_output==1) then
      write(60,'(14E26.16)') time/T_char, x_N/x_0, h_FC/H_char, rho_FC/rho_0, &
                             u_FC/U_char, n_aFC, n_sFC, T_FC/T_0, &
                             Fr_FC, Ri_FC, rho_FC/rho_a, phi_aFC, phi_sFC, x_N/x_0-1.d0
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif

    if(non_dim_output==0) then
      if( dabs(x_N-x_N_old_out_front)/dt_out_front/U_char > steady_criteria ) then
        time_steady = tfinal
      else
        if( time < time_steady) then
          time_steady = time
        endif
      endif
    else
      if( dabs(x_N-x_N_old_out_front)/dt_out_front/U_char/T_char > steady_criteria ) then
        time_steady = tfinal*T_char
      else
        if( time < time_steady) then
          time_steady = time
        endif
      endif
    endif
    x_N_old_out_front = x_N

  endif

  !*** dense (H) & deposit (D) ***
  if(n==0) then
    write(70,*)'# t, x_NH, x_ND, x_NH-x_0, x_ND-x_0'
  endif

  if(non_dim_output==0) then
    write(70,'(5E26.16)') time, x_NH, x_ND, x_NH-x_0, x_ND-x_0
  elseif(non_dim_output==1) then
    write(70,'(5E26.16)') time/T_char, x_NH/x_0, x_ND/x_0, x_NH/x_0-1.d0, x_ND/x_0-1.d0
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif

  if(non_dim_output==0) then
    if( dabs(x_NH-x_NH_old_out_front)/dt_out_front/U_char > steady_criteria ) then
      time_H_steady = tfinal
    else
      if( time < time_H_steady) then
        time_H_steady = time
      endif
    endif
  else
    if( dabs(x_NH-x_NH_old_out_front)/dt_out_front/U_char/T_char > steady_criteria ) then
      time_H_steady = tfinal*T_char
    else
      if( time < time_H_steady) then
        time_H_steady = time
      endif
    endif
  endif
  x_NH_old_out_front = x_NH



  write(*,*)'*** output_front ***'
  write(*,*)'time step n=',n
  if(non_dim_output==0) then
    write(*,*)'time=',time
    write(*,*)'x_N=',x_N
    write(*,*)'x_NH=',x_NH
    write(*,*)'x_ND=',x_ND
  elseif(non_dim_output==1) then
    write(*,*)'time/T_char=',time/T_char
    write(*,*)'x_N/x_0=',x_N/x_0
    write(*,*)'x_NH/x_0=',x_NH/x_0
    write(*,*)'x_ND/x_0=',x_ND/x_0
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif

  write(*,*)'i_front_L=',i_front_L
  write(*,*)'i_front_H=',i_front_H
  write(*,*)'i_front_D=',i_front_D

  if(lift_off_FC==1) then
    write(*,*)'lift_off_FC=',lift_off_FC
  endif
  if(lift_off==1) then
    write(*,*)'i_liftoff=',i_liftoff
    write(*,*)'FC_lift_off=',FC_liftoff
  endif
  write(*,*)''



  return
end subroutine output_front
