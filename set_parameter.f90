module parameter

  implicit none
  double precision, parameter :: pi = 6.d0*dasin(0.5d0)
  double precision, parameter :: eta_g = 1.d-5
  !*** common parameter ***
  integer, protected :: output_type, output_term, output_point, non_dim_output
  integer, protected :: configuration_type, flux_type, Cd_type
  integer, protected :: u0h0_type, x0_type
  integer, protected :: mx, msource
  integer, protected :: i_start
  double precision, protected :: left_boundary_point
  double precision, protected :: x_0, H_x0, h0_x0, y_0, y0_x0
  double precision, protected :: tfinal, dt_out ! [s or -]
  double precision, protected :: dt_out_front, supply_time ! [s or -]
  double precision, protected :: cfl, grav, theta_slope
  double precision, protected :: slope_distance ! [m or -]
  double precision, protected :: x_output_1, x_output_2 ! [m or -]
  double precision, protected :: x_output_3, x_output_4 ! [m or -]
  double precision, protected :: rho_s, rho_a, T_a, R_a, R_v, T_w, rho_w, T_l
  double precision, protected :: a_v
  double precision, protected :: C_s, C_pa, C_pv, C_p0, C_pw
  double precision, protected :: steady_criteria
  double precision, protected :: p_a0, dx
  double precision, protected :: U_char, T_char, H_char
  !*** dilute parameter ***
  integer, protected :: FC_type
  integer, protected :: thermal_type, entrainment_type, ent_mechanic_type
  integer, protected :: Ws_type, geometry_type
  integer, protected :: potential_type, tau_c_type
  double precision, protected :: cfl_FC, Ws_directly, diameter
  double precision, protected :: n_a0, n_s0, phi_a0, phi_s0, rho_0, T_0, n_w0, phi_w0
  double precision, protected :: u_0, h_0, M_dot_0
  double precision, protected :: Ri_0, Fr_0
  double precision, protected :: Fr_N0, C_dc
  double precision, protected :: Re_0, eta_0
  double precision, protected :: p_v0, T_boil0
  !*** dense parameter ***
  integer, protected :: no_dense_layer, D_type
  integer, protected :: basal_resistance_type, pressure_gradient_type
  integer, protected :: particle_momentum_type
  double precision, protected :: D_Ws, D_directly
  double precision, protected :: eps, delta, phi_sH, phi_sD, C_db
  double precision, protected :: rho_H, rho_gH, g_H, mu
  double precision, protected :: Re_H0, eta_H0, h_Hchar
  double precision, protected :: h_H0, u_H0

  namelist /common_list/ output_type,output_term,output_point,non_dim_output,&
                        configuration_type,&
                        flux_type,Cd_type,&
                        u0h0_type,x0_type,mx,msource,&
                        left_boundary_point,x_0,H_x0,h0_x0,&
                        y_0,y0_x0,&
                        tfinal,dt_out,dt_out_front,&
                        supply_time,&
                        cfl,&
                        grav,theta_slope,slope_distance,&
                        x_output_1,x_output_2,x_output_3,x_output_4,&
                        rho_s,p_a0,T_a,&
                        R_a,R_v,C_s,C_pa,C_pv,&
                        rho_w,T_w,C_pw,T_l,&
                        steady_criteria
  namelist /dilute_list/ FC_type,&
                        thermal_type,entrainment_type,ent_mechanic_type,&
                        Ws_type,geometry_type,&
                        potential_type,tau_c_type,&
                        n_s0,n_a0,T_0,&
                        cfl_FC,Fr_N0,&
                        Ws_directly,diameter,&
                        C_dc,&
                        h_0,u_0,&
                        M_dot_0,Ri_0
  namelist /dense_list/ no_dense_layer, D_type, basal_resistance_type, &
                       pressure_gradient_type,particle_momentum_type, &
                       D_Ws, D_directly, &
                       eps,delta,phi_sH,phi_sD, &
                       h_H0,u_H0, &
                       C_db

contains
  subroutine set_parameter()
    open(1,file='PARAM.nml',form='formatted')
    read(1,nml=common_list)
    read(1,nml=dilute_list)
    read(1,nml=dense_list)
    close(1)

    !*** Left-Boundary point ***
    if(left_boundary_point<0.d0 .and. configuration_type==2) then
      write(*,*)'*** ERROR *** (in set_parameter.f90)'
      write(*,*)'  left_boundary_point =', left_boundary_point
      write(*,*)'  configuration_type =', configuration_type
      write(*,*)'Set left_boundary_point>=0 or configuration_type=1'
      stop
    endif

    !*** From degree to rad ***
    theta_slope = theta_slope * pi / 180.d0
    delta = delta * pi / 180.d0
    mu = dtan(delta)

    !*** Ambient conditions ***
    a_v = 0.d0
    rho_a = ((1.d0-a_v)*R_a + a_v*R_v)*T_a/p_a0
    rho_a = 1.d0 / rho_a

    n_w0 = 0.d0
    C_p0 = n_s0*C_s + n_a0*C_pa + (1.d0-n_s0-n_a0)*C_pv
    rho_0 = n_s0/rho_s + ( n_a0*R_a + (1.d0-n_s0-n_a0)*R_v ) * T_0 / p_a0
    rho_0 = 1.d0 / rho_0
    if(rho_0/rho_a<=1.d0) then
      write(*,*)'*** ERROR *** (rho_0<=rho_a)'
      write(*,*)'  rho_0 =', rho_0
      write(*,*)'  rho_a =', rho_a
      stop
    endif
    
    phi_s0 = rho_0*n_s0 / rho_s
    phi_a0 = rho_0*n_a0 * R_a*T_0/p_a0
    phi_w0 = rho_0*n_w0 / rho_w

    if(u0h0_type == 1) then
      Fr_0 = 1.d0 / dsqrt(Ri_0)
      if(x0_type==0) then
        if(configuration_type==1) then
          h_0 = ( M_dot_0 / (y_0*rho_0*Fr_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          H_char = ( M_dot_0 / (y_0*rho_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
        elseif(configuration_type==2) then
          h_0 = ( M_dot_0 / (2.d0*pi*x_0*rho_0*Fr_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          H_char = ( M_dot_0 / (2.d0*pi*x_0*rho_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
        endif
        h_0 = (h_0)**(1.d0/3.d0)
        H_char = (H_char)**(1.d0/3.d0)
        u_0 = Fr_0 * dsqrt((rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) * h_0)
      elseif(x0_type==1) then
        if(configuration_type==1) then
          H_char = ( H_x0*M_dot_0 / (y0_x0*rho_0) )**2.d0
          H_char = H_char / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          H_char = (H_char)**(1.d0/5.d0)
          x_0 = H_char / H_x0
          y_0 = x_0 * y0_x0
          h_0 = ( M_dot_0 / (y_0*rho_0*Fr_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          h_0 = (h_0)**(1.d0/3.d0)
        elseif(configuration_type==2) then
          H_char = ( H_x0*M_dot_0 / (2.d0*pi*rho_0) )**2.d0
          H_char = H_char / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          H_char = (H_char)**(1.d0/5.d0)
          x_0 = H_char / H_x0
          h_0 = ( M_dot_0 / (2.d0*pi*x_0*rho_0*Fr_0) )**2.d0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) )
          h_0 = (h_0)**(1.d0/3.d0)
        endif
        u_0 = Fr_0 * dsqrt((rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) * h_0)
      elseif(x0_type==2) then
        if(configuration_type==1) then
          write(*,*)'*** UNDER COSTRUCTION ***'
          write(*,*)'  x0_type =', x0_type
          write(*,*)'  configuration_type =', configuration_type
          stop
        elseif(configuration_type==2) then
          h_0 = ( h0_x0*M_dot_0 / (2.d0*pi*rho_0) )**2.d0
          h_0 = h_0 / ( (rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) * Fr_0*Fr_0 )
          h_0 = (h_0)**(1.d0/5.d0)
          x_0 = h_0 / h0_x0
          u_0 = Fr_0 * dsqrt((rho_0-rho_a)/rho_0 * grav * dcos(theta_slope) * h_0)
        else
          write(*,*)'*** ERROR ***'
          write(*,*)'  x0_type =', x0_type
          write(*,*)'  configuration_type =', configuration_type
          stop
        endif
      else
        write(*,*)'*** ERROR ***'
        write(*,*)'   x0_type =', x0_type
        stop
      endif
    elseif(u0h0_type==2 .and. x0_type==0 .and. configuration_type==1) then
      H_char = h_0
      M_dot_0 = y_0*rho_0*u_0*h_0
      Ri_0 = (rho_0-rho_a)*grav*dcos(theta_slope)*h_0 / (rho_0*u_0*u_0)
      Fr_0 = 1.d0 / dsqrt(Ri_0)
    else
      write(*,*)'*** ERROR (in set_parameter.f90) ***'
      write(*,*)'  u0h0_type =', u0h0_type
      write(*,*)'  x0_type =', x0_type
      write(*,*)'  configuration_type =', configuration_type
      stop
    endif


    rho_gH = p_a0 / (R_a*T_0)
    rho_H = phi_sH*rho_s + (1.d0-phi_sH)*rho_gH
    g_H = grav*(rho_H-rho_a)/rho_H


    if(x0_type==2) then
      H_char = h_0
      U_char = u_0
    else
      U_char = dsqrt(H_char*grav*(rho_0-rho_a)/rho_0)
    endif
    T_char = x_0 / U_char

    if(non_dim_output==1) then ! from non-dim to dim for numerical simulation
      supply_time = supply_time * T_char
      left_boundary_point = left_boundary_point * x_0
      x_output_1 = x_output_1 * x_0
      x_output_2 = x_output_2 * x_0
      x_output_3 = x_output_3 * x_0
      x_output_4 = x_output_4 * x_0
    endif

    eta_0 = eta_g * (1.d0+2.5d0*phi_s0+10.05d0*phi_s0**2+0.00273d0*dexp(16.6d0*phi_s0))
    eta_H0 = eta_g * (1.d0+2.5d0*phi_sH+10.05d0*phi_sH**2+0.00273d0*dexp(16.6d0*phi_sH))
    h_Hchar = (n_s0*rho_0)/(phi_sH*rho_s)*Ws_func(rho_0,T_0,n_s0,n_a0,n_w0)*dcos(theta_slope)
    h_Hchar = h_Hchar - phi_sD/phi_sH*D_func(rho_H,T_0,phi_sH*rho_s/rho_H,0.d0,0.d0)*dcos(theta_slope)
    h_Hchar = h_Hchar * T_char
    if(h_Hchar<h_H0) h_Hchar = h_H0
    Re_0 = rho_0 * U_char * H_char / eta_0
    Re_H0 = rho_H * U_char * h_Hchar / eta_H0
    if(Cd_type==1) then
      C_dc = 25.d-3 * Re_0**(-0.2d0)
      C_db = 25.d-3 * Re_H0**(-0.2d0)
    endif
    if(h_Hchar<=0.d0) then
      C_db = 0.d0
      Re_H0 = 0.d0
    endif
    h_Hchar = (n_s0*rho_0)/(phi_sH*rho_s)*Ws_func(rho_0,T_0,n_s0,n_a0,n_w0)*dcos(theta_slope)*T_char
    if(h_Hchar<h_H0) h_Hchar = h_H0

    dx = x_0 / dble(msource) ! msource[grid]:x_0[m or -] = 1[grid]:dx[m or -]
    if(non_dim_output==1) then
      slope_distance = slope_distance * x_0 ! from non-dim to dim
    endif

    i_start = msource + 1

    if(Ws_type==3.and.Ws_directly<0.d0) then
      write(*,*)'*** ERROR *** (in set_parameter.f90)'
      write(*,*)'  Ws_type =', Ws_type
      write(*,*)'  Ws_directly =', Ws_directly
      stop
    elseif(Ws_type==1.or.Ws_type==2) then
      if(diameter<=0.d0) then
        write(*,*)'*** ERROR *** (in set_parameter.f90)'
        write(*,*)'  Ws_type =', Ws_type
        write(*,*)'  diameter =', diameter
        stop
      endif
    endif

    T_boil0 = 0.d0

    return
  end subroutine set_parameter



  double precision function Ent_func(Fr)

    implicit none

    double precision, intent(in) :: Fr
    double precision :: Ri, E

    if(Fr<=0.d0) then
      Ri = 0.d0
    else
      Ri = 1.d0 / Fr / Fr
    endif

    if(entrainment_type==0) then
      E = 0.d0
    elseif(entrainment_type==1) then ! Parker (1987)
    if(Fr<=0.d0) then
      E = 0.d0
    else
      E = 0.075d0 / dsqrt( 1.d0 + 718.d0*Ri**(2.4d0) )
    endif
    elseif(entrainment_type==2) then ! Turner (1986)
      E = (0.08d0*Fr*Fr - 0.1d0) / (Fr*Fr+5.d0)
      E = max(E, 0.d0)
    elseif(entrainment_type==3) then ! Johnson & Hogg (2013)
      if(Fr<=0.d0) then
        E = 0.d0
      else
        E = 0.075d0 / ( 1.d0 + 27.d0*Ri )
      endif
    else
      write(*,*)'*** ERROR *** (entrainment_type=???)'
      write(*,*)'   entrainment_type=', entrainment_type
      stop
    endif

    ent_func = E

    return
  end function ent_func




  double precision function Ws_func(rho,T,n_s,n_a,n_w)

    implicit none

    double precision, intent(in) :: rho, T, n_s, n_a, n_w
    double precision :: rho_g, R_g
    double precision :: Ws

    if(Ws_type==1) then ! Bursik & Woods (1996)
      Ws = dsqrt( rho_s * grav * diameter / rho )
    elseif(Ws_type==2) then ! Doyle et al. (2008;2010;2011)
      R_g = ( n_a*R_a + (1.d0-n_s-n_w-n_a)*R_v ) / (1.d0-n_s-n_w)
      rho_g = p_a0 / (R_g*T)
      Ws = dsqrt(4.d0*(rho_s-rho_g)*grav*diameter / (3.d0*rho_g))
    elseif(Ws_type==3) then ! Ws is directly given.
      Ws = Ws_directly
    else
      write(*,*)'*** ERROR *** (Ws_type=???)'
      stop
    endif

    Ws_func = Ws


    return
  end function Ws_func




  double precision function D_func(rho,T,n_s,n_a,n_w)

    implicit none

    double precision, intent(in) :: rho, T, n_s, n_a, n_w
    double precision :: D, Ws

    if(D_type==0) then ! D/Ws is directly given.
      Ws = Ws_func(rho,T,n_s,n_a,n_w)
!      if(Ws==0.d0) then
!        D = 0.d0
!      else
        D = D_Ws * Ws
!      endif
    elseif(D_type==1) then ! D is directly given.
      D = D_directly
    else
      write(*,*)'*** ERROR *** (D_type=???)'
      stop
    endif

    D_func = D


    return
  end function D_func




end module parameter
