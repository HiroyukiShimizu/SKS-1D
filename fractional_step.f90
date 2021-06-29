subroutine fractional_step(n,dt,time)

  use variable
  use parameter
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: time
  double precision, intent(inout) :: dt

  integer, parameter :: niter = 10
  integer :: i, it, i_front_max
  double precision :: s1, s2, s3, s4, s5, s6, s5_work, s6_work ! source term
  double precision :: tau_c, tau_b, pressure_gradient, sigma, sigma_b
  double precision :: Ws, Ent, D
  double precision :: f1, f2, f3, f4, f5, f6
  double precision :: A1, A2, A3, A4, A5, A6
  double precision :: h_FC, rho_FC, u_FC
  double precision :: na_FC, ns_FC, phia_FC, phis_FC, nw_FC, phiw_FC
  double precision :: T_FC, enthal_FC, Cp_FC, aaa_FC, Fr_FC
  double precision :: Tboil_FC
  double precision :: nv_work, pv_work, Tboil_work, T_work, Cp_work
  double precision :: ns_work, na_work, nw_work, rhoh_work, enthal_work

  do it = 1, niter

    !**************
    !*  DELETE
    !**************
    !*** dilute (L) ***
    Q_star(:,:) = 0.d0
    F(:,:) = 0.d0
    h_star(:) = 0.d0
    rho_star(:) = 0.d0
    u_star(:) = 0.d0
    na_star(:) = 0.d0
    ns_star(:) = 0.d0
    nw_star(:) = 0.d0
    phia_star(:) = 0.d0
    phis_star(:) = 0.d0
    phiw_star(:) = 0.d0
    aaa_star(:) = 0.d0
    Fr_star(:) = 0.d0
    T_star(:) = 0.d0
    T_boil_star(:) = 0.d0
    enthal_star(:) = 0.d0
    Cp_star(:) = 0.d0
    c(:) = 0.d0
    !*** dense (H) ***
    QH_star(:,:) = 0.d0
    FH(:,:) = 0.d0
    hH_star(:) = 0.d0
    uH_star(:) = 0.d0
    aaaH_star(:) = 0.d0
    FrH_star(:) = 0.d0
    interact(:,:) = 0.d0
    !*** deposit (D) ***
    dhD(:) = 0.d0
    dhD_from_H(:) = 0.d0


    !*************************************
    !***  1st step
    !*************************************

    !*** left boundary condition (solid wall) ***
    z_b(i_start-1) = z_b(i_start)
    z_c(i_start-1) = z_c(i_start)
    if(time+dt <= supply_time) then
      rho(i_start-1) = rho_0
      h(i_start-1) = h_0
      u(i_start-1) = u_0
    else
      rho(i_start-1) = rho(i_start)
      h(i_start-1) = h(i_start)
      u(i_start-1) = -u(i_start)
    endif

    !*** dilute (L) ***
    if(lift_off==0) then
      do i = i_start, i_front_L
        if(rho(i)*0.d0/=0.d0) then
          write(*,*)'i =', i
          write(*,*)'rho(i) =', rho(i)
          stop
        endif
        s1 = 0.d0
        s2 = 0.d0
        s3 = 0.d0
        s4 = 0.d0
        s5 = 0.d0
        s6 = 0.d0

        !*** slope ***
        if(h(i)>H_char*eps.and.rho(i)>rho_a) then
          s4 = s4 + (rho(i)-rho_a)*grav*h(i)*dsin(theta(i))
        endif

        !*** basal resistance ***
        if(h(i)>H_char*eps.and.rho(i)>rho_a) then
          if(tau_c_type==1) then
            tau_c = C_dc * rho(i) * (u(i)-uH(i)) * dabs(u(i)-uH(i))
          elseif(tau_c_type==2) then
            tau_c = C_dc * rho(i) * u(i) * dabs(u(i))
          else
            write(*,*)'*** ERROR *** (tau_c_type)'
            write(*,*)'  tau_c_type =', tau_c_type
            stop
          endif
          if(tau_c*0.d0/=0.d0) then
            write(*,*)'i =',i
            write(*,*)'FC =',FC
            write(*,*)'tau_c =',tau_c
            write(*,*)'rho(i) =',rho(i)
            write(*,*)'u(i) =',u(i)
            write(*,*)'uH(i) =',uH(i)
            tau_c = 0.d0
          endif
          s4 = s4 - tau_c
          if(i==FC) then
            interact(i,2) = interact(i,2) - tau_c/rho_H * dx_FC/dx
          else
            interact(i,2) = interact(i,2) - tau_c/rho_H
          endif
        endif

        !*** particle settling ***
        if(h(i)>H_char*eps.and.rho(i)>rho_a) then
          Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
          s1 = s1 - rho(i)*Ws*dcos(theta(i))*ns(i)
          s3 = s3 - rho(i)*Ws*dcos(theta(i))*ns(i)
          s4 = s4 - rho(i)*Ws*dcos(theta(i))*ns(i) * u(i)
          s5 = s5 - rho(i)*Ws*dcos(theta(i))*ns(i) * C_s*T(i)
          if(potential_type==0) then
            s5 = s5 + rho(i)*Ws*dcos(theta(i))*ns(i) * 0.5d0*grav*h(i)*dcos(theta(i))
          endif
          if(i==FC) then
            interact(i,1) = interact(i,1) &
                            + rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH/rho_s * dx_FC/dx
            if(particle_momentum_type==0) then
              interact(i,2) = interact(i,2) &
                              + rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH*u(i)/rho_s*dx_FC/dx
            endif
          else
            interact(i,1) = interact(i,1) &
                            + rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH/rho_s
            if(particle_momentum_type==0) then
              interact(i,2) = interact(i,2) &
                              + rho(i)*Ws*dcos(theta(i))*ns(i)/phi_sH*u(i)/rho_s
            endif
          endif
        endif

        !*** geometry ***
        if(geometry_type==1.and.h(i)>H_char*eps.and.rho(i)>rho_a) then
          s4 = s4 - (rho(i)-rho_a)*grav*h(i)*dcos(theta(i))*(z_c(i)-z_c(i-1))/dx
        endif

        !*** entrainment ***
        if(entrainment_type/=0.and.h(i)>H_char*eps.and.rho(i)>rho_a) then
          Ent = Ent_func(Fr(i))
          s1 = s1 + rho_a*dabs(u(i))*Ent
          s2 = s2 + (1.d0-a_v)*rho_a*dabs(u(i))*Ent
          if(ent_mechanic_type==1) then
            s5 = s5 + rho_a*dabs(u(i))*Ent * ((1.d0-a_v)*C_pa+a_v*C_pv) * T_a
          elseif(ent_mechanic_type==0) then
            s5 = s5 + rho_a*dabs(u(i))*Ent * ( ((1.d0-a_v)*C_pa+a_v*C_pv)*T_a + 0.5d0*u(i)*u(i) )
          else
            write(*,*)'*** ERROR (in fractional_step.f90) ***'
            write(*,*)'   ent_mechanic_type =', ent_mechanic_type
            stop
          endif
          if(potential_type==0) then
            s5 = s5 + rho_a*dabs(u(i))*Ent * 0.5d0*grav*h(i)*dcos(theta(i))
          endif
        endif

        !*** No considering potential energy ***
        if(potential_type==1) then
          s5 = s5 + rho(i)*grav*u(i)*h(i)*dcos(theta(i))*(0.5d0*(h(i)-h(i-1))+z_c(i)-z_c(i-1))/dx
        elseif(potential_type==3) then
          s5 = s5 + 0.5d0*u(i)*grav*dcos(theta(i))*((rho(i)-rho_a)*h(i)*h(i)-(rho(i-1)-rho_a)*h(i-1)*h(i-1))/dx
        endif

        !*** axisymmetric ***  
        if(configuration_type==2) then
          s1 = s1 - Q(i,1)*u(i)/x(i)
          s2 = s2 - Q(i,2)*u(i)/x(i)
          s3 = s3 - Q(i,3)*u(i)/x(i)
          s4 = s4 - Q(i,4)*u(i)/x(i)
          s5 = s5 - Q(i,5)*u(i)/x(i)
          s6 = s6 - Q(i,6)*u(i)/x(i)
        endif

        !*** interaction of pressure gradient ***
        if(pressure_gradient_type==1.and.h(i)>H_char*eps.and.rho(i)>rho_a) then
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
          interact(i,2) = interact(i,2) - pressure_gradient
        endif

        if(interact(i,1)*0.d0/=0.d0 .or. interact(i,2)*0.d0/=0.d0) then
          write(*,*)'i =',i
          write(*,*)'FC =',FC
          write(*,*)'interact(i,1) =',interact(i,1)
          write(*,*)'interact(i,2) =',interact(i,2)
          write(*,*)'pressure_gradient =',pressure_gradient
          write(*,*)'tau_c/rho_H =',tau_c/rho_H
          write(*,*)'rho(i) =',rho(i)
          stop
        endif

        !**************************************************
        Q_star(i,1) = Q(i,1) + dt*s1
        Q_star(i,2) = Q(i,2) + dt*s2
        Q_star(i,3) = Q(i,3) + dt*s3
        Q_star(i,4) = Q(i,4) + dt*s4
        if(thermal_type/=0) Q_star(i,5) = Q(i,5) + dt*s5
        Q_star(i,6) = Q(i,6) + dt*s6
        !**************************************************

      enddo

      T(i_start-1) = T(i_start)
      enthal(i_start-1) = enthal(i_start)
      call cons_to_prim(Q_star,h_star,rho_star,u_star,&
                        na_star,ns_star,nw_star,phia_star,phis_star,phiw_star,&
                        aaa_star,Fr_star,T_star,enthal_star,Cp_star,&
                        u,T,enthal,dt,i_front_L,1,theta,T_boil_star)

      !*** LIFT OFF ??? ***
      do i = i_front_L, i_start, -1

        if(rho_star(i)*0.d0/=0.d0) then
          write(*,*)'i =',i
          write(*,*)'rho_star(i) =',rho_star(i)
          stop
        endif

        if(rho_star(i)<=rho_a*(1.d0+1.d-4)) then

          if(R_v==R_a.and.thermal_type==0) then
            write(*,*)'i=', i
            write(*,*)'rho_star(i)=',rho_star(i)
            stop
          else
            if(i+1==FC .and. dx_FC==0.d0) then
              Q_star(i,:) = Q_star(FC,:)
              FC = i
            else
              !*** LIFTOFF region ***
              h_star(i) = eps*H_char
              na_star(i) = 1.d0 - a_v
              ns_star(i) = 0.d0
              nw_star(i) = 0.d0
              T_star(i) = T_a
              T_boil_star(i) = 0.d0
              rho_star(i) = rho_a
              aaa_star(i) = 0.d0
              phia_star(i) = na_star(i)*rho_a / ( p_a0/(R_a*T_a) )
              phis_star(i) = 0.d0
              phiw_star(i) = 0.d0
              Cp_star(i) = C_pa
              enthal_star(i) = C_pa * T_a
              Q_star(i,1) = rho_star(i)*h_star(i)
              Q_star(i,2) = na_star(i)*rho_star(i)*h_star(i)
              Q_star(i,3) = 0.d0
              Q_star(i,4) = 0.d0
              Q_star(i,5) = rho_star(i)*enthal_star(i)*h_star(i)
              Q_star(i,6) = 0.d0
              u_star(i) = 0.d0
              Fr_star(i) = 0.d0
              if(i==FC) then
                dx_FC = 0.d0
                lift_off_FC = 1
              endif
            endif
          endif

        endif

      enddo

    endif

    !*** dense (H) ***
    if(no_dense_layer==1) then
      do i = i_start, i_front_L
        dhD(i) = dt*interact(i,1)
      enddo
    else
      i_front_max = max(i_front_L,i_front_H)
      do i = i_start, i_front_max
        s1 = 0.d0
        s2 = 0.d0

        !*** interaction ***
        s1 = s1 + interact(i,1)
        s2 = s2 + interact(i,2)

        !*** sedimentation ***
        D = D_func(rho_H,T_0,phi_sH*rho_s/rho_H,0.d0,0.d0)
        s1 = s1 - D*dcos(theta(i)) * phi_sD/phi_sH
        s2 = s2 - D*uH(i)*dcos(theta(i)) * phi_sD/phi_sH

        !*** axisymmetric ***
        if(configuration_type==2) then
          s1 = s1 - hH(i)*uH(i)/x(i)
          s2 = s2 - hH(i)*uH(i)*uH(i)/x(i)
        endif

        !*************************
        !*** judgement of mass ***
        !*************************
        if(hH(i)+dt*s1<=h_Hchar*eps) then
          dhD(i) = ( hH(i) + dt*interact(i,1) ) * phi_sH/phi_sD
          dhD_from_H(i) = hH(i) * phi_sH/phi_sD
          hH_star(i) = 0.d0
          uH_star(i) = 0.d0
          aaaH_star(i) = 0.d0
          FrH_star(i) = 0.d0
        else
          dhD(i) = dt*D*dcos(theta(i))
          dhD_from_H(i) = dt*D*dcos(theta(i))
          hH_star(i) = hH(i) + dt*s1
          if(hH_star(i)*0.d0/=0.d0) then
            write(*,*)'i =', i
            write(*,*)'hH(i) =', hH(i)
            write(*,*)'hH_star(i) =', hH_star(i)
            write(*,*)'s1 =', s1
          endif

          !*** slope ***
          s2 = s2 + g_H*hH(i)*dsin(theta(i))

          !*** geometry ***
          s2 = s2 - g_H*hH(i)*dcos(theta(i)) * (z_b(i)-z_b(i-1))/dx

          !*** basal resistance ***
          if(basal_resistance_type==0) then
            tau_b = 0.d0
          elseif(basal_resistance_type==1) then
            tau_b = mu*(rho_H-rho_a)*grav*hH(i)*dcos(theta(i))*dsign(1.d0,uH(i))
          elseif(basal_resistance_type==2) then
            tau_b = rho_H*C_db*uH(i)*dabs(uH(i))
          else
            write(*,*)'*** ERROR *** (basal_resistance_type=???)'
            write(*,*)'  basal_resistance_type=',basal_resistance_type
            stop
          endif
          !*****************************
          !*** judgement of momentum ***
          !*****************************
          sigma_b = dabs(dt*tau_b/rho_H)
          sigma = dabs(uH(i)*hH(i) + dt*s2)
          if(sigma_b>=sigma) then
            uH_star(i) = 0.d0
            aaaH_star(i) = 0.d0
            FrH_star(i) = 0.d0
          else
            s2 = s2 - tau_b/rho_H
            uH_star(i) = uH(i)*hH(i) + dt*s2
            uH_star(i) = uH_star(i)/hH_star(i)
            aaaH_star(i) = dsqrt(g_H*hH_star(i)*dcos(theta(i)))
            FrH_star(i) = dabs(uH_star(i)) / aaaH_star(i)
          endif

        endif

      enddo

      !*** Artificial-Bed model ***
      do i = i_start, mx
        if(hH_star(i)>0.d0) then
          QH_star(i,1) = hH_star(i) + h_Hchar * eps
          QH_star(i,2) = QH_star(i,1) * uH_star(i)
        elseif(hH_star(i)<=0.d0) then
          QH_star(i,1) = h_Hchar * eps
          QH_star(i,2) = 0.d0
        else
          write(*,*)'*** ERROR *** (hH_star(i)=???)'
          write(*,*)'   i=',i
          write(*,*)'   hH_star(i)=',hH_star(i)
          write(*,*)'   hH_star(i-1)=',hH_star(i-1)
          write(*,*)'   hH_star(i+1)=',hH_star(i+1)
          write(*,*)'   interact(i,1)=',interact(i,1)
          stop
        endif
      enddo

    endif


    !*************************************
    !***  2nd step
    !*************************************
    !*** dilute (L) ***
    if(lift_off==0) then
      if(time+dt <= supply_time) then
        h_star(i_start-1) = h_0
        rho_star(i_start-1) = rho_0
        u_star(i_start-1) = u_0
        na_star(i_start-1) = n_a0
        ns_star(i_start-1) = n_s0
        nw_star(i_start-1) = n_w0
        phia_star(i_start-1) = phi_a0
        phis_star(i_start-1) = phi_s0
        phiw_star(i_start-1) = phi_w0
        aaa_star(i_start-1) = dsqrt(h_0*grav*dcos(theta(i_start-1))*(rho_0-rho_a)/rho_0)
        T_star(i_start-1) = T_0
        T_boil_star(i_start-1) = T_boil0
        enthal_star(i_start-1) = C_p0 * T_0
      else
        h_star(i_start-1) = h_star(i_start)
        rho_star(i_start-1) = rho_star(i_start)
        u_star(i_start-1) = - u_star(i_start)
        na_star(i_start-1) = na_star(i_start)
        ns_star(i_start-1) = ns_star(i_start)
        nw_star(i_start-1) = nw_star(i_start)
        phia_star(i_start-1) = phia_star(i_start)
        phis_star(i_start-1) = phis_star(i_start)
        phiw_star(i_start-1) = phiw_star(i_start)
        aaa_star(i_start-1) = aaa_star(i_start)
        T_star(i_start-1) = T_star(i_start)
        T_boil_star(i_start-1) = T_boil_star(i_start)
        enthal_star(i_start-1) = enthal_star(i_start)
      endif

      if(dx_FC<=0.d0) then
      !*** numerical flux ***
        if(flux_type==1) then
          call hll_flux_L(FC-1, F, h_star, rho_star, u_star, &
                          na_star, ns_star, nw_star, aaa_star, enthal_star, theta)
        elseif(flux_type==2) then
          call hllc_flux_L(FC-1, F, h_star, rho_star, u_star, &
                           na_star, ns_star, nw_star, aaa_star, enthal_star, theta)
        else
          write(*,*)'*** ERROR *** (flux_type)'
          write(*,*)'   flux_type=',flux_type
          stop
        endif

        call godunov_flux_FC(h_star(FC-1),rho_star(FC-1),u_star(FC-1),&
                             na_star(FC-1),ns_star(FC-1),nw_star(FC-1),aaa_star(FC-1),&
                             Fr_star(FC-1),enthal_star(FC-1),f1,f2,f3,f4,f5,f6,&
                             theta(FC-1))

        F(FC,1) = f1
        F(FC,2) = f2
        F(FC,3) = f3
        F(FC,4) = f4
        F(FC,5) = f5
        F(FC,6) = f6

      else
      !*** numerical flux ***
        if(flux_type==1) then
          call hll_flux_L(FC, F, h_star, rho_star, u_star, &
                          na_star, ns_star, nw_star, aaa_star, enthal_star, theta)
        elseif(flux_type==2) then
          call hllc_flux_L(FC, F, h_star, rho_star, u_star, &
                           na_star, ns_star, nw_star, aaa_star, enthal_star, theta)
        else
          write(*,*)'*** ERROR *** (flux_type)'
          write(*,*)'   flux_type=',flux_type
          stop
        endif
        f1 = F(FC,1)
        f2 = F(FC,2)
        f3 = F(FC,3)
        f4 = F(FC,4)
        f5 = F(FC,5)
        f6 = F(FC,6)

      endif

      A1 = dx_FC*Q_star(FC,1) + dt*f1
      A2 = dx_FC*Q_star(FC,2) + dt*f2
      A3 = dx_FC*Q_star(FC,3) + dt*f3
      A4 = dx_FC*Q_star(FC,4) + dt*f4
      A5 = dx_FC*Q_star(FC,5) + dt*f5
      A6 = dx_FC*Q_star(FC,6) + dt*f6

      if(A4>0.d0) then
        call cons_to_prim_FC(A1, A2, A3, A4, A5, A6, dt, &
                             h_FC, rho_FC, u_FC, &
                             na_FC, ns_FC, nw_FC, phia_FC, phis_FC, phiw_FC, &
                             T_FC, enthal_FC, Cp_FC, aaa_FC, Fr_FC, &
                             theta(FC), Tboil_FC)
      endif
      if(lift_off_FC==1) then
        !*** LIFTOFF region ***
        dx_FC_new = 0.d0
        h_FC = eps*H_char
        na_FC = 1.d0 - a_v
        ns_FC = 0.d0
        nw_FC = 0.d0
        T_FC = T_a
        Tboil_FC = 0.d0
        rho_FC = rho_a
        aaa_FC = 0.d0
        phia_FC = na_FC*rho_a / ( p_a0/(R_a*T_a) )
        phis_FC = 0.d0
        phiw_FC = 0.d0
        Cp_FC = C_pa
        enthal_FC = C_pa * T_a
        u_FC = 0.d0
        Fr_FC = 0.d0
      endif

      if(dx_FC_new-dx_FC>cfl*dx) then
        if(it==1) then
          dt = cfl_FC*dx / ( (dx_FC_new-dx_FC)/dt )
        else
          dt = 0.9d0*dt
        endif
      else
        exit
      endif

    endif

    if(lift_off==1) exit

    if(it>=niter) then
      write(*,*)'*** ERROR *** (dx_FC vs dt)'
      write(*,*)'  it=', it
      write(*,*)'  niter=', niter
      stop
    endif

  enddo

  if(lift_off==0) then

    do i = i_start, FC-1
      !**************************************************
      Q(i,1) = Q_star(i,1) - dt/dx * ( F(i+1,1) - F(i,1) )
      Q(i,2) = Q_star(i,2) - dt/dx * ( F(i+1,2) - F(i,2) )
      Q(i,3) = Q_star(i,3) - dt/dx * ( F(i+1,3) - F(i,3) )
      Q(i,4) = Q_star(i,4) - dt/dx * ( F(i+1,4) - F(i,4) )
      if(thermal_type/=0) Q(i,5) = Q_star(i,5) - dt/dx * ( F(i+1,5) - F(i,5) )
      Q(i,6) = Q_star(i,6) - dt/dx * ( F(i+1,6) - F(i,6) )
      !**************************************************
    enddo

    if(time+dt <= supply_time) then
      T_star(i_start-1) = T_0
      enthal_star(i_start-1) = C_p0*T_0
    else
      T_star(i_start-1) = T_star(i_start)
      enthal_star(i_start-1) = enthal_star(i_start)
    endif

    call cons_to_prim(Q,h,rho,u,na,ns,nw,phia,phis,phiw,aaa,Fr,T,enthal,Cp,&
                      u_star,T_star,enthal_star,dt,FC-1,2,theta,T_boil)

    h(FC) = h_FC
    rho(FC) = rho_FC
    u(FC) = u_FC
    na(FC) = na_FC
    ns(FC) = ns_FC
    nw(FC) = nw_FC
    phia(FC) = phia_FC
    phis(FC) = phis_FC
    phiw(FC) = phiw_FC
    aaa(FC) = aaa_FC
    if(aaa(FC)*0.d0/=0.d0) then
      write(*,*)'*** ERROR *** (aaa(FC)=NaN)'
      write(*,*)'   aaa(FC)=',aaa(FC)
      stop
    endif
    Fr(FC) = Fr_FC
    T(FC) = T_FC
    T_boil(FC) = Tboil_FC
    enthal(FC) = enthal_FC
    Cp(FC) = Cp_FC
    Q(FC,1) = rho_FC * h_FC
    Q(FC,2) = na_FC * rho_FC * h_FC
    Q(FC,3) = ns_FC * rho_FC * h_FC
    Q(FC,4) = rho_FC * u_FC * h_FC
    if(thermal_type/=0) Q(FC,5) = rho_FC * enthal_FC * h_FC
    Q(FC,6) = nw_FC * rho_FC * h_FC

    !*** LIFT OFF ??? ***
    do i = FC-1, i_start, -1

      if(rho(i)*0.d0/=0.d0) then
        write(*,*)'i =', i
        write(*,*)'rho(i) =', rho(i)
        write(*,*)'rho_a =', rho_a
        stop
      endif

      if(rho(i)<=rho_a*(1.d0+1.d-4)) then
        if(R_v==R_a.and.thermal_type==0) then
          write(*,*)'i=',i
          write(*,*)'rho(i)=',rho(i)
          stop
        else
          if(i+1==FC .and. dx_FC_new==0.d0) then
            Q(i,:) = Q(FC,:)
            FC = i
            write(*,*)'i =', i
            write(*,*)'FC =', FC
          else
            if(i==FC-1) then
              write(*,*)'i =', i
              write(*,*)'FC =', FC
            endif

            !*** LIFTOFF region ***
            h(i) = eps*H_char
            na(i) = 1.d0 - a_v
            ns(i) = 0.d0
            nw(i) = 0.d0
            T(i) = T_a
            T_boil(i) = 0.d0
            rho(i) = rho_a
            aaa(i) = 0.d0
            phia(i) = na(i)*rho_a / ( p_a0/(R_a*T_a) )
            phis(i) = 0.d0
            phiw(i) = 0.d0
            Cp(i) = C_pa
            enthal(i) = C_pa * T_a
            Q(i,1) = rho(i)*h(i)
            Q(i,2) = na(i)*rho(i)*h(i)
            Q(i,3) = 0.d0
            Q(i,4) = 0.d0
            Q(i,5) = rho(i)*enthal(i)*h(i)
            Q(i,6) = 0.d0
            u(i) = 0.d0
            Fr(i) = 0.d0
          endif
        endif
      endif

    enddo

  endif


  !*** dense (H) ***
  if(no_dense_layer==0) then
    !*** left boundary condition ***
    if(h_H0>0.d0.and.time+dt<=supply_time) then
      QH_star(i_start-1,1) = h_H0 + h_Hchar*eps
      QH_star(i_start-1,2) = u_H0 * (h_H0 + h_Hchar*eps)
    else
      QH_star(i_start-1,1) = QH_star(i_start,1)
      QH_star(i_start-1,2) = - QH_star(i_start,2)
    endif
    !*** right boundary condition (solid wall) ***
    QH_star(mx+1,1) = QH_star(mx,1)
    QH_star(mx+1,2) = - QH_star(mx,2)

    !*** numerical flux ***
    if(flux_type==1) then
      call hll_flux_H(FH, QH_star,theta)
    elseif(flux_type==2) then
      call hllc_flux_H(FH, QH_star,theta)
    else
      write(*,*)'*** ERROR *** (flux_type)'
      write(*,*)'   flux_type=',flux_type
      stop
    endif

    do i = i_start, mx
      !**************************************************
      QH(i,1) = QH_star(i,1) - dt/dx * ( FH(i+1,1) - FH(i,1) )
      QH(i,2) = QH_star(i,2) - dt/dx * ( FH(i+1,2) - FH(i,2) )
      !**************************************************

      !*** Artificial-Bed model ***
      if(QH(i,1)<=h_Hchar*eps) then
        QH(i,1) = h_Hchar*eps
        QH(i,2) = 0.d0
        hH(i) = 0.d0
        uH(i) = 0.d0
        aaaH(i) = 0.d0
        FrH(i) = 0.d0
      else
        hH(i) = QH(i,1)-h_Hchar*eps
        uH(i) = QH(i,2)/QH(i,1)
        aaaH(i) = dsqrt(g_H*hH(i)*dcos(theta(i)))
        FrH(i) = dabs(uH(i)) / aaaH(i)
      endif
    enddo

  endif


  return
end subroutine fractional_step

