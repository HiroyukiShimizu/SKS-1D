subroutine output(iframe,n,time,dt)

  use parameter
  use variable
  implicit none

  integer, intent(in) :: iframe, n
  double precision, intent(in) :: time, dt

  integer :: i, nstp, ipos, idigit
  double precision :: z_f, Ri, xxx, RiH, p, E, Eu, nv, Sa
  double precision :: inertia, hydrostatic, tauc_drag, supplied_particle
  double precision :: pressure_gradient, geometry, taub_drag, sedimentation
  double precision :: axisymmetric, slope, inertia_r
  double precision :: tau_b, Ws, tau_c, D
  double precision :: massflux, massflux_r
  double precision :: entrainment, Ent
  double precision :: phi_sH_from_L
  double precision :: ent_mechanic, particle, gas, Cp_effect
  double precision :: inertia_nondim, entrainment_nondim, ent_mechanic_nondim
  double precision :: particle_nondim, gas_nondim, Cp_effect_nondim
  double precision :: massflux_nondim, supplied_particle_nondim
  double precision :: sedimentation_nondim, axisymmetric_nondim, massflux_r_nondim
  double precision :: hydrostatic_nondim, tauc_drag_nondim
  double precision :: pressure_gradient_nondim
  double precision :: geometry_nondim, taub_drag_nondim
  double precision :: slope_nondim, inertia_r_nondim
  double precision :: no_potential, no_potential_nondim
  double precision :: ent_potential, ent_potential_nondim
  double precision :: particle_potential, particle_potential_nondim
  double precision :: condensation, condensation_nondim
  character(len=10) :: fname1, fname2, fname3, fname4, fname5
  character(len=10) :: fname10, fname20, fname30, fname40, fname50


  phi_sH_from_L = 0.d0
  
  !*** dilute (L) ***
  fname1 = 'prim.Lxxxx'
  nstp = iframe
  do ipos=10,7,-1
    idigit = mod(nstp,10)
    fname1(ipos:ipos) = char(ichar('0')+idigit)
    nstp = nstp / 10
  enddo
  open(unit=10,file=fname1,status='unknown',form='formatted')
  if(non_dim_output==0) then
    write(10,10) time
  elseif(non_dim_output==1) then
    write(10,10) time/T_char
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif
  10 format('### t= ',f7.2)
  write(10,*)'# x h rho u na ns T Fr Ri rho/rho_a z_f enthal p/p_a0'
  write(10,*)'# phia phis phi_sH_from_L E Eu T/T_a x-x_0 T_boil 1-ns-nw ns/ns0 nw phiw nv'
  if(lift_off==0) then
    do i = i_start, i_front_L+1
      z_f = z_c(i) + h(i)
      if(Fr(i)<=0.d0) then
        Ri = 0.d0
        E = 0.d0
        Eu = 0.d0
      elseif(Fr(i)>0.d0) then
        Ri = 1.d0 / Fr(i) / Fr(i)
        E = Ent_func(Fr(i))
        Eu = E*u(i)
      endif
      if(i<i_front_L) then
        xxx = x(i)
      elseif(i>=i_front_L) then
        xxx = x_N
      endif
      if(i<=i_front_L) then
        p = T(i)*(na(i)*R_a+(1.d0-ns(i)-na(i)-nw(i))*R_v) / (1.d0/rho(i) - ns(i)/rho_s - nw(i)/rho_w)
        nv = 1.d0-ns(i)-nw(i)-na(i)
      else
        p = 0.d0
        nv = 0.d0
      endif

      if(non_dim_output==0) then
        write(10,'(26E26.16)') xxx,h(i),rho(i),u(i),na(i),ns(i),T(i),Fr(i),Ri,&
                               rho(i)/rho_a,z_f,enthal(i),p/p_a0,&
                               phia(i),phis(i),phi_sH_from_L,E,Eu,T(i)/T_a,xxx-x_0,&
                               T_boil(i),1.d0-ns(i)-nw(i),ns(i)/n_s0,nw(i),phiw(i),&
                               nv
      elseif(non_dim_output==1) then
        write(10,'(26E26.16)') xxx/x_0,h(i)/H_char,rho(i)/rho_0,u(i)/U_char,&
                               na(i),ns(i),T(i)/T_0,Fr(i),Ri,&
                               rho(i)/rho_a,z_f/H_char,enthal(i)/(C_p0*T_0),p/p_a0,&
                               phia(i),phis(i),phi_sH_from_L,E,Eu/U_char,T(i)/T_a,xxx/x_0-1.d0,&
                               T_boil(i)/T_0,1.d0-ns(i)-nw(i),ns(i)/n_s0,nw(i),phiw(i),&
                               nv
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    enddo
  endif
  write(10,*)' '
  close(10)

  fname2 = 'cons.Lxxxx'
  nstp = iframe
  do ipos=10,7,-1
    idigit = mod(nstp,10)
    fname2(ipos:ipos) = char(ichar('0')+idigit)
    nstp = nstp / 10
  enddo
  open(unit=20,file=fname2,status='unknown',form='formatted')
  if(non_dim_output==0) then
    write(20,20) time
  elseif(non_dim_output==1) then
    write(20,20) time/T_char
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif
  20 format('### t= ',f7.2)
  write(20,*)'# x rho*h na*rho*h ns*rho*h rho*u*h rho*enthal*h rho*h*u na*rho*h*u ns*rho*h*u rho*u*h*u rho*enthal*h*u'
  write(20,*)'# ns*rho*h*u*u (2*pi*)ns*rho*h*u*r ns*rho*h*u*u*r x-x_0 [ns*rho*h*u*r/(ns0*rho0*h0*u0*r0)]'
  write(20,*)'# nw*rho*h nw*rho*h*u'
  if(lift_off==0) then
    do i = i_start, i_front_L+1
      if(i<i_front_L) xxx = x(i)
      if(i>=i_front_L) xxx = x_N
      if(non_dim_output==0) then
        write(20,'(18E26.16)') xxx,Q(i,1),Q(i,2),Q(i,3),Q(i,4),Q(i,5),&
                              Q(i,1)*u(i),Q(i,2)*u(i),Q(i,3)*u(i),Q(i,4)*u(i),Q(i,5)*u(i),&
                              Q(i,3)*u(i)*u(i),2.d0*pi*Q(i,3)*u(i)*x(i),Q(i,3)*u(i)*u(i)*x(i),xxx-x_0, &
                              Q(i,3)*u(i)*x(i)/(n_s0*rho_0*U_char*H_char*x_0), &
                              Q(i,6),Q(i,6)*u(i)
      elseif(non_dim_output==1) then
        write(20,'(17E26.16)') xxx/x_0, Q(i,1)/(rho_0*H_char), Q(i,2)/(rho_0*H_char), &
                              Q(i,3)/(rho_0*H_char), Q(i,4)/(rho_0*U_char*H_char), &
                              Q(i,5)/(rho_0*C_p0*T_0*H_char), &
                              Q(i,1)*u(i)/(rho_0*H_char*U_char), Q(i,2)*u(i)/(rho_0*H_char*U_char), &
                              Q(i,3)*u(i)/(rho_0*H_char*U_char), Q(i,4)*u(i)/(rho_0*U_char*H_char*U_char), &
                              Q(i,5)*u(i)/(rho_0*C_p0*T_0*H_char*U_char),&
                              Q(i,3)*u(i)*u(i)/(rho_0*U_char*H_char*U_char),&
                              Q(i,3)*u(i)*x(i)/(rho_0*U_char*H_char*x_0),&
                              Q(i,3)*u(i)*u(i)*x(i)/(rho_0*U_char*H_char*U_char*x_0),xxx/x_0-1.d0, &
                              Q(i,6)/(rho_0*H_char), Q(i,6)/(rho_0*U_char*H_char)
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    enddo
  endif
  write(20,*)' '
  close(20)


  if(output_term==1) then

    fname40 = 'moti.Lxxxx'
    nstp = iframe
    do ipos=10,7,-1
      idigit = mod(nstp,10)
      fname40(ipos:ipos) = char(ichar('0')+idigit)
      nstp = nstp / 10
    enddo
    open(unit=140,file=fname40,status='unknown',form='formatted')
    if(non_dim_output==0) then
      write(140,140) time
    elseif(non_dim_output==1) then
      write(140,140) time/T_char
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif
    140 format('### t= ',f7.2)
    write(140,*)'# x inertia hydrostatic slope geometry tauc_drag entrainment x-x_0'
    if(lift_off==0) then
      z_c(i_start-1) = z_c(i_start)
      uH(i_start-1) = -uH(i_start)
      do i = i_start, i_front_L+1
        if(i<i_front_L) xxx = x(i)
        if(i>=i_front_L) xxx = x_N
        if(i<=i_front_L) then
          inertia = u(i)*(u(i)-u(i-1))/dx
          hydrostatic = - 0.5d0*grav*dcos(theta(i)) / (rho(i)*h(i)) 
          hydrostatic = hydrostatic * ((rho(i)-rho_a)*h(i)*h(i) - (rho(i-1)-rho_a)*h(i-1)*h(i-1)) / dx
          slope = grav*dsin(theta(i))*(rho(i)-rho_a)/rho(i)
          geometry = - grav*dcos(theta(i)) * (rho(i)-rho_a)/rho(i) * (z_c(i)-z_c(i-1))/dx
          if(tau_c_type==1) then
            tau_c = C_dc * rho(i) * (u(i)-uH(i)) * dabs(u(i)-uH(i))
          elseif(tau_c_type==2) then
            tau_c = C_dc * rho(i) * u(i) * dabs(u(i))
          else
            tau_c = 0.d0
          endif
          tauc_drag = - tau_c / (rho(i)*h(i))
          Ent = Ent_func(Fr(i))
          entrainment = -rho_a*dabs(u(i))*Ent*u(i) / (rho(i)*h(i))
        else
          inertia = 0.d0 
          hydrostatic = 0.d0
          slope = 0.d0
          geometry = 0.d0
          tauc_drag = 0.d0
          entrainment = 0.d0
        endif
        if(non_dim_output==0) then
          write(140,'(8E26.16)') xxx,inertia,hydrostatic,slope,geometry,tauc_drag,entrainment,xxx-x_0
        elseif(non_dim_output==1) then
          inertia_nondim = inertia * x_0/(U_char*U_char)
          hydrostatic_nondim = hydrostatic * x_0/(U_char*U_char)
          slope_nondim = slope * x_0/(U_char*U_char)
          geometry_nondim = geometry * x_0/(U_char*U_char)
          tauc_drag_nondim = tauc_drag * x_0/(U_char*U_char)
          entrainment_nondim = entrainment * x_0/(U_char*U_char)
          write(140,'(8E26.16)') xxx/x_0,inertia_nondim,hydrostatic_nondim,slope_nondim,&
                                 geometry_nondim,tauc_drag_nondim,entrainment_nondim,xxx/x_0-1.d0
        else
          write(*,*)'*** ERROR *** (non_dim_output=???)'
          stop
        endif
      enddo
    endif
    write(140,*)' '
    close(140)

    fname50 = 'temp.Lxxxx'
    nstp = iframe
    do ipos=10,7,-1
      idigit = mod(nstp,10)
      fname50(ipos:ipos) = char(ichar('0')+idigit)
      nstp = nstp / 10
    enddo
    open(unit=150,file=fname50,status='unknown',form='formatted')
    if(non_dim_output==0) then
      write(150,150) time
    elseif(non_dim_output==1) then
      write(150,150) time/T_char
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif
    150 format('### t= ',f7.2)
    write(150,*)'# x inertia entrainment particle gas ent_mechanic Cp_effect no_potential x-x_0  condensation'
    if(lift_off==0) then
      z_c(i_start-1) = z_c(i_start)
      Cp(i_start-1) = C_p0
      h(i_start-1) = h_0
      rho(i_start-1) = rho_0
      do i = i_start, i_front_L+1
        if(i<i_front_L) xxx = x(i)
        if(i>=i_front_L) xxx = x_N
        if(i<=i_front_L) then
          inertia = u(i)*(T(i)-T(i-1))/dx
          Ent = Ent_func(Fr(i))
          entrainment = rho_a*dabs(u(i))*Ent
          entrainment = entrainment * (((1.d0-a_v)*C_pa+a_v*C_pv)*T_a-enthal(i))
          entrainment = entrainment / (Cp(i)*rho(i)*h(i))
          if(ent_mechanic_type==1) then
            ent_mechanic = 0.d0
          elseif(ent_mechanic_type==0) then
            ent_mechanic = rho_a*dabs(u(i))*Ent * 0.5d0*u(i)*u(i) / (Cp(i)*rho(i)*h(i))
          endif
          if(h(i)>H_char*eps.and.rho(i)>rho_a.and.i<=i_front_L) then
            Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
            particle = -rho(i)*Ws*dcos(theta(i))*ns(i) * (C_s*T(i)-enthal(i)) / (Cp(i)*rho(i)*h(i))
            gas = 0.d0
            if(potential_type==0) then
              particle_potential = rho(i)*Ws*dcos(theta(i))*ns(i) * 0.5d0*grav*h(i)*dcos(theta(i)) / (Cp(i)*rho(i)*h(i))
            else
              particle_potential = 0.d0
            endif
            condensation = c(i) * (2.5d6-(C_pw-C_pv)*(T(i)-273.d0)) / (Cp(i)*rho(i)*h(i))
          else
            particle = 0.d0
            gas = 0.d0
            particle_potential = 0.d0
          endif
          if(potential_type==0) then
            ent_potential = rho_a*dabs(u(i))*Ent * 0.5d0*grav*h(i)*dcos(theta(i)) / (Cp(i)*rho(i)*h(i))
            no_potential = 0.d0
          elseif(potential_type==1) then
            ent_potential = 0.d0
            no_potential = rho(i)*grav*u(i)*h(i)*dcos(theta(i))*(0.5d0*(h(i)-h(i-1))+z_c(i)-z_c(i-1))/dx / (Cp(i)*rho(i)*h(i))
          elseif(potential_type==3) then
            ent_potential = 0.d0
            no_potential = 0.5d0*u(i)*grav*dcos(theta(i))*((rho(i)-rho_a)*h(i)*h(i)-(rho(i-1)-rho_a)*h(i-1)*h(i-1))/dx
            no_potential = no_potential/(Cp(i)*rho(i)*h(i))
          else
            ent_potential = 0.d0
            no_potential = 0.d0
          endif
          if(thermal_type==2) then
            if(n==0) then
              Cp_effect = - T(i)/Cp(i) * ( u(i)*(Cp(i)-Cp(i-1))/dx )
            else
              Cp_effect = - T(i)/Cp(i) * ( (Cp(i)-Cp_old(i))/dt +u(i)*(Cp(i)-Cp(i-1))/dx )
            endif
          else
            Cp_effect = 0.d0
          endif
        else
          inertia = 0.d0 
          entrainment = 0.d0
          ent_mechanic = 0.d0
          particle = 0.d0
          gas = 0.d0
          Cp_effect = 0.d0
          no_potential = 0.d0
          ent_potential = 0.d0
          particle_potential = 0.d0
          condensation = 0.d0
        endif
        if(non_dim_output==0) then
          write(150,'(12E26.16)') xxx,inertia,entrainment,particle,gas,ent_mechanic,Cp_effect,&
                                 no_potential,ent_potential,particle_potential,xxx-x_0,condensation
        elseif(non_dim_output==1) then
          inertia_nondim = inertia * x_0/(T_0*U_char)
          entrainment_nondim = entrainment * x_0/(T_0*U_char)
          particle_nondim = particle * x_0/(T_0*U_char)
          gas_nondim = gas * x_0/(T_0*U_char)
          ent_mechanic_nondim = ent_mechanic * x_0/(T_0*U_char)
          Cp_effect_nondim = Cp_effect * x_0/(T_0*U_char)
          no_potential_nondim = no_potential * x_0/(T_0*U_char)
          ent_potential_nondim = ent_potential * x_0/(T_0*U_char)
          particle_potential_nondim = particle_potential * x_0/(T_0*U_char)
          condensation = condensation * x_0/(T_0*U_char)
          write(150,'(12E26.16)') xxx/x_0,inertia_nondim,entrainment_nondim,particle_nondim,&
                                 gas_nondim,ent_mechanic_nondim,Cp_effect_nondim,&
                                 no_potential_nondim,ent_potential_nondim,particle_potential_nondim,&
                                 xxx/x_0-1.d0,condensation
        else
          write(*,*)'*** ERROR *** (non_dim_output=???)'
          stop
        endif
      enddo
    endif
    write(150,*)' '
    close(150)

  endif


  if(no_dense_layer==0) then
    !*** dense (H) ***
    fname3 = 'prim.Hxxxx'
    nstp = iframe
    do ipos=10,7,-1
      idigit = mod(nstp,10)
      fname3(ipos:ipos) = char(ichar('0')+idigit)
      nstp = nstp / 10
    enddo
    open(unit=30,file=fname3,status='unknown',form='formatted')
    if(non_dim_output==0) then
      write(30,30) time
    elseif(non_dim_output==1) then
      write(30,30) time/T_char
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif
    30 format('### t= ',f7.2)
    write(30,*)'# x  hH  uH  FrH  RiH  z_c  Sa  x-x_0'
    do i = i_start, i_front_H+1
      if(FrH(i)<=0.d0) RiH = 0.d0
      if(FrH(i)>0.d0) RiH = 1.d0 / FrH(i) / FrH(i)
      if(Ws_type==1.or.Ws_type==2) then
        Sa = rho_s * (diameter*uH(i)/hH(i))**2 / ( (rho_s-rho_a)*grav*hH(i)*dcos(theta(i)) )
      else
        Sa = 0.d0
      endif
      if(non_dim_output==0) then
        write(30,'(8E26.16)') x(i),hH(i),uH(i),FrH(i),RiH,z_c(i),Sa,x(i)-x_0
      elseif(non_dim_output==1) then
        write(30,'(8E26.16)') x(i)/x_0,hH(i)/H_char,uH(i)/U_char,FrH(i),RiH,z_c(i)/H_char,Sa,x(i)/x_0-1.d0
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    enddo
    write(30,*)' '
    close(30)

    fname4 = 'cons.Hxxxx'
    nstp = iframe
    do ipos=10,7,-1
      idigit = mod(nstp,10)
      fname4(ipos:ipos) = char(ichar('0')+idigit)
      nstp = nstp / 10
    enddo
    open(unit=40,file=fname4,status='unknown',form='formatted')
    if(non_dim_output==0) then
      write(40,40) time
    elseif(non_dim_output==1) then
      write(40,40) time/T_char
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif
    40 format('### t= ',f7.2)
    write(40,*)'# x  hH  uH*hH  uH*hH*uH  uH*hH/dx  uH*hH*uH/dx'
    write(40,*)'# uH*hH*phi_sH*rho_s  uH*uH*hH*phi_sH*rho_s  uH*hH*phi_sH*rho_s*r  uH*uH*hH*phi_sH*rho_s*r'
    write(40,*)'# x-x_0'
    do i = i_start, i_front_H+1
      if(non_dim_output==0) then
        write(40,'(11E26.16)') x(i),QH(i,1),QH(i,2),QH(i,2)*uH(i),QH(i,2)/dx,QH(i,2)*uH(i)/dx, &
                              QH(i,2)*phi_sH*rho_s, QH(i,2)*uH(i)*phi_sH*rho_s, &
                              QH(i,2)*phi_sH*rho_s*x(i), QH(i,2)*uH(i)*phi_sH*rho_s*x(i), x(i)-x_0
      elseif(non_dim_output==1) then
        write(40,'(11E26.16)') x(i)/x_0,QH(i,1)/H_char,QH(i,2)/(U_char*H_char), &
                              QH(i,2)*uH(i)/(U_char*H_char*U_char), &
                              QH(i,2)/dx * x_0/(U_char*H_char), &
                              QH(i,2)*uH(i)/dx * x_0/(U_char*H_char*U_char), &
                              QH(i,2)*phi_sH*rho_s/(U_char*H_char*rho_0), &
                              QH(i,2)*uH(i)*phi_sH*rho_s/(U_char*U_char*H_char*rho_0), &
                              QH(i,2)*phi_sH*rho_s*x(i)/(U_char*H_char*rho_0*x_0), &
                              QH(i,2)*uH(i)*phi_sH*rho_s*x(i)/(U_char*U_char*H_char*rho_0*x_0), &
                              x(i)/x_0-1.d0
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    enddo
    write(40,*)' '
    close(40)


    if(output_term==1) then

      fname10 = 'mome.Hxxxx'
      fname20 = 'moti.Hxxxx'
      fname30 = 'mass.Hxxxx'
      nstp = iframe
      do ipos=10,7,-1
        idigit = mod(nstp,10)
        fname10(ipos:ipos) = char(ichar('0')+idigit)
        fname20(ipos:ipos) = char(ichar('0')+idigit)
        fname30(ipos:ipos) = char(ichar('0')+idigit)
        nstp = nstp / 10
      enddo
      open(unit=110,file=fname10,status='unknown',form='formatted')
      open(unit=120,file=fname20,status='unknown',form='formatted')
      open(unit=130,file=fname30,status='unknown',form='formatted')
      if(non_dim_output==0) then
        write(110,110) time
        write(120,120) time
        write(130,130) time
      elseif(non_dim_output==1) then
        write(110,110) time/T_char
        write(120,120) time/T_char
        write(130,130) time/T_char
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
      110 format('### t= ',f7.2)
      write(110,*)'# x inertia hydrostaic tauc_drag supplied_particle pressure_grad'
      write(110,*)'# geometry taub_drag sedimen axi slope inertia_r x-x_0'
      120 format('### t= ',f7.2)
      write(120,*)'# x inertia hydrostaic tauc_drag supplied_particle pressure_grad geometry taub_drag slope x-x_0'
      130 format('### t= ',f7.2)
      write(130,*)'# x massflux supplied_particle sedimen axi massflux_r x-x_0'
      hH(i_start-1) = hH(i_start)
      uH(i_start-1) = -uH(i_start)
      z_b(i_start-1) = z_b(i_start)
      D = D_func(rho_H,T_0,phi_sH*rho_s/rho_H,0.d0,0.d0)
      do i = i_start, i_front_H+1
        massflux = (uH(i)*hH(i) - uH(i-1)*hH(i-1)) / dx
        if(h(i)>H_char*eps.and.rho(i)>rho_a.and.i<=i_front_L) then
          Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
          if(i==FC) then
            supplied_particle = rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH/rho_s*dx_FC/dx
          else
            supplied_particle = rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH/rho_s
          endif
        else
          supplied_particle = 0.d0
        endif
        if(configuration_type==1) then
          axisymmetric = 0.d0
        elseif(configuration_type==2) then
          axisymmetric = - hH(i)*uH(i)/x(i)
        endif
        massflux_r = massflux - axisymmetric
        if(non_dim_output==0) then
          write(130,'(7E26.16)') x(i),massflux,supplied_particle,-D*dcos(theta(i))*phi_sD/phi_sH,axisymmetric,massflux_r,x(i)-x_0
        elseif(non_dim_output==1) then
          massflux_nondim = massflux * x_0/(U_char*H_char)
          supplied_particle_nondim = supplied_particle * x_0/(U_char*H_char)
          sedimentation_nondim = - D*dcos(theta(i)) * phi_sD/phi_sH * x_0/(U_char*H_char)
          axisymmetric_nondim = axisymmetric * x_0/(U_char*H_char)
          massflux_r_nondim = massflux_r * x_0/(U_char*H_char)
          write(130,'(7E26.16)') x(i)/x_0,massflux_nondim,supplied_particle_nondim, &
                                 sedimentation_nondim,axisymmetric_nondim,massflux_r_nondim, &
                                 x(i)/x_0-1.d0
        else
          write(*,*)'*** ERROR *** (non_dim_output=???)'
          stop
        endif

        inertia = (uH(i)*uH(i)*hH(i) - uH(i-1)*uH(i-1)*hH(i-1)) / dx
        hydrostatic = - 0.5d0*g_H*dcos(theta(i)) * (hH(i)*hH(i) - hH(i-1)*hH(i-1)) / dx
        if(h(i)>H_char*eps.and.rho(i)>rho_a) then
          if(tau_c_type==1) then
            tau_c = C_dc * rho(i) * (u(i)-uH(i)) * dabs(u(i)-uH(i))
          elseif(tau_c_type==2) then
            tau_c = C_dc * rho(i) * u(i) * dabs(u(i))
          else
            tau_c = 0.d0
          endif
          if(i==FC) then
            tauc_drag = - tau_c/rho_H * dx_FC/dx
          else
            tauc_drag = - tau_c/rho_H
          endif
        else
          tauc_drag = 0.d0
        endif
        supplied_particle = supplied_particle*u(i)
        if(pressure_gradient_type==1.and.h(i)>H_char*eps.and.rho(i)>rho_a.and.i<=i_front_L) then
          if(i==FC) then
            pressure_gradient = hH(i)*grav*dcos(theta(i))/rho_H/dx &
                                *( - (rho(i-1)-rho_a)*h(i-1) )
          else
            pressure_gradient = hH(i)*grav*dcos(theta(i))/rho_H/dx &
                                *((rho(i)-rho_a)*h(i) - (rho(i-1)-rho_a)*h(i-1))
          endif
          if(pressure_gradient*0.d0/=0.d0) then
            write(*,*)'i =',i
            write(*,*)'FC =',FC
            write(*,*)'pressure_gradient =',pressure_gradient
            pressure_gradient = 0.d0
          endif
          pressure_gradient = - pressure_gradient
        else
          pressure_gradient = 0.d0
        endif
        geometry = - g_H*hH(i)*dcos(theta(i)) * (z_b(i)-z_b(i-1))/dx
        if(basal_resistance_type==0) then
          tau_b = 0.d0
        elseif(basal_resistance_type==1) then
          tau_b = mu*(rho_H-rho_a)*grav*hH(i)*dcos(theta(i))*dsign(1.d0,uH(i))
        elseif(basal_resistance_type==2) then
          tau_b = rho_H*C_db*uH(i)*dabs(uH(i))
        endif
        taub_drag = - tau_b / rho_H
        sedimentation = - D*uH(i)*dcos(theta(i)) * phi_sD/phi_sH
        slope = g_H*hH(i)*dsin(theta(i))
        if(configuration_type==1) then
          axisymmetric = 0.d0
        elseif(configuration_type==2) then
          axisymmetric = - hH(i)*uH(i)*uH(i)/x(i)
        endif
        inertia_r = inertia - axisymmetric
        if(non_dim_output==0) then
          write(110,'(13E26.16)') x(i),inertia,hydrostatic,tauc_drag,supplied_particle,pressure_gradient,&
                                  geometry,taub_drag,sedimentation,axisymmetric,slope,inertia_r,x(i)-x_0
        elseif(non_dim_output==1) then
          inertia_nondim = inertia * x_0/(U_char*U_char*H_char)
          hydrostatic_nondim = hydrostatic * x_0/(U_char*U_char*H_char)
          tauc_drag_nondim = tauc_drag * x_0/(U_char*U_char*H_char)
          supplied_particle_nondim = supplied_particle * x_0/(U_char*U_char*H_char)
          pressure_gradient_nondim = pressure_gradient * x_0/(U_char*U_char*H_char)
          geometry_nondim = geometry * x_0/(U_char*U_char*H_char)
          taub_drag_nondim = taub_drag * x_0/(U_char*U_char*H_char)
          sedimentation_nondim = sedimentation * x_0/(U_char*U_char*H_char)
          axisymmetric_nondim = axisymmetric * x_0/(U_char*U_char*H_char)
          slope_nondim = slope * x_0/(U_char*U_char*H_char)
          inertia_r_nondim = inertia_r * x_0/(U_char*U_char*H_char)
          write(110,'(13E26.16)') x(i)/x_0,inertia_nondim,hydrostatic_nondim,tauc_drag_nondim,&
                                  supplied_particle_nondim,pressure_gradient_nondim,&
                                  geometry_nondim,taub_drag_nondim,sedimentation_nondim,&
                                  axisymmetric_nondim,slope_nondim,inertia_r_nondim,&
                                  x(i)/x_0-1.d0
        else
          write(*,*)'*** ERROR *** (non_dim_output=???)'
          stop
        endif

        if(i<i_front_H+1) then
          inertia = uH(i) * ( uH(i) - uH(i-1) ) / dx
          hydrostatic = hydrostatic / hH(i)
          tauc_drag = tauc_drag / hH(i)
          if(u(i)/=0.d0) then
            supplied_particle = supplied_particle * (u(i)-uH(i))/u(i)/hH(i)
          else
            supplied_particle = 0.d0
          endif
          pressure_gradient = pressure_gradient / hH(i)
          geometry = geometry / hH(i)
          taub_drag = taub_drag / hH(i)
          slope = slope / hH(i)
        else
          inertia = uH(i-1) * ( uH(i) - uH(i-1) ) / dx
          hydrostatic = hydrostatic / QH(i,1)
          tauc_drag = tauc_drag / QH(i,1)
          if(u(i)/=0.d0) then
            supplied_particle = supplied_particle * (u(i)-uH(i))/u(i)/QH(i,1)
          else
            supplied_particle = 0.d0
          endif
          pressure_gradient = pressure_gradient / QH(i,1)
          geometry = geometry / QH(i,1)
          taub_drag = taub_drag / QH(i,1)
          slope = slope / QH(i,1)
        endif
        if(non_dim_output==0) then
          write(120,'(10E26.16)') x(i),inertia,hydrostatic,tauc_drag,supplied_particle,pressure_gradient,&
                                  geometry,taub_drag,slope,x(i)-x_0
        elseif(non_dim_output==1) then
          inertia_nondim = inertia * x_0/(U_char*U_char)
          hydrostatic_nondim = hydrostatic * x_0/(U_char*U_char)
          tauc_drag_nondim = tauc_drag * x_0/(U_char*U_char)
          supplied_particle_nondim = supplied_particle * x_0/(U_char*U_char)
          pressure_gradient_nondim = pressure_gradient * x_0/(U_char*U_char)
          geometry_nondim = geometry * x_0/(U_char*U_char)
          taub_drag_nondim = taub_drag * x_0/(U_char*U_char)
          slope_nondim = slope * x_0/(U_char*U_char)
          write(120,'(10E26.16)') x(i)/x_0,inertia_nondim,hydrostatic_nondim,tauc_drag_nondim,&
                                 supplied_particle_nondim,pressure_gradient_nondim,&
                                 geometry_nondim,taub_drag_nondim,slope_nondim,x(i)/x_0-1.d0
        else
          write(*,*)'*** ERROR *** (non_dim_output=???)'
          stop
        endif
      enddo
      write(110,*)' '
      close(110)
      write(120,*)' '
      close(120)
      write(130,*)' '
      close(130)

    endif

  endif



  !*** deposit (D) ***
  fname5 = 'prim.Dxxxx'
  nstp = iframe
  do ipos=10,7,-1
    idigit = mod(nstp,10)
    fname5(ipos:ipos) = char(ichar('0')+idigit)
    nstp = nstp / 10
  enddo
  open(unit=50,file=fname5,status='unknown',form='formatted')
  if(non_dim_output==0) then
    write(50,50) time
  elseif(non_dim_output==1) then
    write(50,50) time/T_char
  else
    write(*,*)'*** ERROR *** (non_dim_output=???)'
    stop
  endif
  50 format('### t= ',f7.2)
  write(50,*)'# x  hD  hD(from H)  x-x_0  phi_sD*rho_s*hD  z_b'
  do i = i_start, i_front_D+1
    if(non_dim_output==0) then
      write(50,'(6E26.16)') x(i),hD(i),hD_from_H(i),x(i)-x_0,phi_sD*rho_s*hD(i),z_b(i)
    elseif(non_dim_output==1) then
      write(50,'(6E26.16)') x(i)/x_0,hD(i)/H_char,hD_from_H(i)/H_char,x(i)/x_0-1.d0,phi_sD*rho_s*hD(i)/(rho_0*H_char),z_b(i)/H_char
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif
  enddo
  write(50,*)' '
  close(50)

  !*** Initial topography z_b0 ***
  if(time==0.d0) then
    write(250,*)'# x  z_b0'
    do i = 1, mx+1
      if(non_dim_output==0) then
        write(250,'(2E26.16)') x(i),z_b0(i)
      elseif(non_dim_output==1) then
        write(250,'(2E26.16)') x(i)/x_0,z_b0(i)/H_char
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    enddo
    write(250,*)' '
    close(250)
  endif

  write(*,*)'*** output ***'
  write(*,*)'time step n=',n
  if(non_dim_output==0) write(*,*)'time=',time
  if(non_dim_output==1) write(*,*)'time/T_char=',time/T_char
  write(*,*)''


  return
end subroutine output
