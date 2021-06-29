program main

  use variable
  use parameter
  implicit none

  integer :: i, n, iframe
  integer :: i_output_1,i_output_2,i_output_3,i_output_4
  double precision :: time, dt, time_out, time_out_front

  open(unit=60,file='front_L.dat',status='unknown',form='formatted')
  open(unit=70,file='front_H.dat',status='unknown',form='formatted')
  open(unit=200,file='mass.dat',status='unknown',form='formatted')
  open(unit=250,file='initial_zb.dat',status='unknown',form='formatted')

  call set_parameter()

  if(output_point==1) then
    open(unit=600,file='point_1_L.dat',status='unknown',form='formatted')
    open(unit=650,file='point_1_H.dat',status='unknown',form='formatted')
    open(unit=700,file='point_2_L.dat',status='unknown',form='formatted')
    open(unit=750,file='point_2_H.dat',status='unknown',form='formatted')
    open(unit=800,file='point_3_L.dat',status='unknown',form='formatted')
    open(unit=850,file='point_3_H.dat',status='unknown',form='formatted')
    open(unit=900,file='point_4_L.dat',status='unknown',form='formatted')
    open(unit=950,file='point_4_H.dat',status='unknown',form='formatted')
  endif

  !*** common ***
  allocate(x(0:mx+1))
  allocate(theta(0:mx+1))

  !*** dilute (L) ***
  allocate(Q(0:mx+1,1:6), Q_star(0:mx+1,1:6)) !--> 1:rho*h 2:na*rho*h 3:ns*rho*h 4:rho*u*h
  allocate(F(0:mx+1,1:6))                     !    5:rho*Cp*T*h 6:nw*rho*h
  allocate(h(0:mx+1), h_star(0:mx+1))
  allocate(rho(0:mx+1), rho_star(0:mx+1))
  allocate(u(0:mx+1), u_star(0:mx+1))
  allocate(na(0:mx+1), na_star(0:mx+1))
  allocate(ns(0:mx+1), ns_star(0:mx+1))
  allocate(nw(0:mx+1), nw_star(0:mx+1))
  allocate(phia(0:mx+1), phia_star(0:mx+1))
  allocate(phis(0:mx+1), phis_star(0:mx+1))
  allocate(phiw(0:mx+1), phiw_star(0:mx+1))
  allocate(aaa(0:mx+1), aaa_star(0:mx+1))
  allocate(Fr(0:mx+1), Fr_star(0:mx+1))
  allocate(T(0:mx+1), T_star(0:mx+1))
  allocate(enthal(0:mx+1), enthal_star(0:mx+1))
  allocate(Cp(0:mx+1), Cp_star(0:mx+1), Cp_old(0:mx+1))
  allocate(z_c(0:mx+1))
  allocate(T_boil(0:mx+1), T_boil_star(0:mx+1))
  allocate(c(0:mx+1))

  !*** dense (H) ***
  allocate(QH(0:mx+1,1:2), QH_star(0:mx+1,1:2)) !--> 1:hH 2:u*hH
  allocate(FH(0:mx+1,1:2))
  allocate(hH(0:mx+1), hH_star(0:mx+1))
  allocate(uH(0:mx+1), uH_star(0:mx+1))
  allocate(aaaH(0:mx+1), aaaH_star(0:mx+1))
  allocate(FrH(0:mx+1), FrH_star(0:mx+1))
  allocate(interact(0:mx+1,1:2))

  !*** deposit (D) ***
  allocate(z_b(0:mx+1),z_b0(0:mx+1))
  allocate(hD(0:mx+1), dhD(0:mx+1))
  allocate(hD_from_H(0:mx+1), dhD_from_H(0:mx+1))

  i_output_1 = 0
  i_output_2 = 0
  i_output_3 = 0
  i_output_4 = 0

  do i = 0, mx+1
    x(i) = dble(i)*dx + left_boundary_point
    if(x(i)>=x_output_1 .and. x(i-1)<=x_output_1) i_output_1 = i ! Point of Profile 2 of PELE
    if(x(i)>=x_output_2 .and. x(i-1)<=x_output_2) i_output_2 = i ! Point of Profile 3 of PELE
    if(x(i)>=x_output_3 .and. x(i-1)<=x_output_3) i_output_3 = i ! Point of Profile 4 of PELE
    if(x(i)>=x_output_4 .and. x(i-1)<=x_output_4) i_output_4 = i ! Point of Profile 5 of PELE

    if(x(i)<=slope_distance) then
      theta(i) = theta_slope
    else
      theta(i) = 0.d0
    endif
  enddo

  call set_initial()

  do i = 0, mx+1
    if(uH(i)*0.d0/=0.d0) then
      write(*,*)'*** After set_initial ***'
      write(*,*)'i =',i
      write(*,*)'FC =',FC
      write(*,*)'uH(i) =',uH(i)
      stop
    endif
  enddo

  iframe = 0
  n = 0
  time = 0.d0
  time_out = 0.d0
  time_out_front = 0.d0
  i_liftoff = mx
  if(non_dim_output==0) then
    time_steady = tfinal
    time_H_steady = tfinal
  else
    time_steady = tfinal*T_char
    time_H_steady = tfinal*T_char
  endif
  x_N_old_out_front = x_N
  x_NH_old_out_front = x_NH

  call output(iframe,n,time,0.d0)
  time_out = time_out + dt_out
  iframe = iframe + 1
  call output_front(n,time)
  call output_mass(n,time)
  if(output_point==1) then
    call output_points(n,time,i_output_1,i_output_2,i_output_3,i_output_4)
  endif
  time_out_front = time_out_front + dt_out_front

  do i = 0, mx+1
    if(uH(i)*0.d0/=0.d0) then
      write(*,*)'*** After output ***'
      write(*,*)'i =',i
      write(*,*)'FC =',FC
      write(*,*)'uH(i) =',uH(i)
      stop
    endif
  enddo

  !********************************************************
  !***   MAIN ROOP
  !********************************************************
  do

    do i = 0, mx+1
      if(uH(i)*0.d0/=0.d0) then
        write(*,*)'*** Before checking cfl_condition ***'
        write(*,*)'i =',i
        write(*,*)'FC =',FC
        write(*,*)'uH(i) =',uH(i)
        uH(i) = 0.d0
      endif
    enddo

    call cfl_condition(dt)

    do i = 0, mx+1
      if(uH(i)*0.d0/=0.d0) then
        write(*,*)'*** After checking cfl_condition ***'
        write(*,*)'i =',i
        write(*,*)'FC =',FC
        write(*,*)'uH(i) =',uH(i)
        stop
      endif
    enddo

    Cp_old(:) = Cp(:)
    call Fractional_step(n,dt,time)


    !*******************
    !*** update
    !*******************
    n = n + 1
    time = time + dt
    x_N_old = x_N
    x_NH_old = x_NH
    if(lift_off==0) then
      !*** dilute (L) ***
      if(FC<=i_start) then

        lift_off = 1
        i_front_L = 0
        x_N = 0.d0

      else

        x_N = dble(FC-1)*dx + dx_FC_new + left_boundary_point
        if(dx_FC_new >= dx) then
          Q(FC+1,1) = Q(FC,1)
          Q(FC+1,2) = Q(FC,2)
          Q(FC+1,3) = Q(FC,3)
          Q(FC+1,4) = Q(FC,4)
          Q(FC+1,5) = Q(FC,5)
          Q(FC+1,6) = Q(FC,6)
          h(FC+1) = h(FC)
          rho(FC+1) = rho(FC)
          u(FC+1) = u(FC)
          na(FC+1) = na(FC)
          ns(FC+1) = ns(FC)
          nw(FC+1) = nw(FC)
          phia(FC+1) = phia(FC)
          phis(FC+1) = phis(FC)
          phiw(FC+1) = phiw(FC)
          aaa(FC+1) = aaa(FC)
          Fr(FC+1) = Fr(FC)
          T(FC+1) = T(FC)
          T_boil(FC+1) = T_boil(FC)
          enthal(FC+1) = enthal(FC)
          Cp(FC+1) = Cp(FC)
          FC = int(x_N/dx) + 1
          dx_FC = dx_FC_new - dx
          if(dx_FC==0.d0) then
            i_front_L = FC-1
          else
            i_front_L = FC
          endif
        else
          if(FC /= int(x_N/dx) + 1) then
            if(dx_FC_new>0.d0) then
              write(*,*)'*** ERROR **** (FC /= int(x_N/dx)+1)'
              write(*,*)'  FC =', FC
              write(*,*)'  int(x_N/dx)+1 =', int(x_N/dx)+1
              write(*,*)'  dx_FC =', dx_FC
              write(*,*)'  dx_FC_new =', dx_FC_new
              stop
            else

            endif
          endif
          dx_FC = dx_FC_new
          if(dx_FC<=0.d0) then
            i_front_L = FC-1
          else
            i_front_L = FC
          endif
        endif

        if(x_N > x_N_max) then
          x_N_max = x_N
        endif
        if(x_N==x_N_old) then
          x_N_steady = x_N
          i_front_steady = int(x_N_steady/dx)
          h_N_steady = h(i_front_steady)
          rho_N_steady = rho(i_front_steady)
          u_N_steady = u(i_front_steady)
          na_N_steady = na(i_front_steady)
          ns_N_steady = ns(i_front_steady)
          nw_N_steady = nw(i_front_steady)
          T_N_steady = T(i_front_steady)
          Fr_N_steady = Fr(i_front_steady)
          enthal_N_steady = enthal(i_front_steady)
          phia_N_steady = phia(i_front_steady)
          phis_N_steady = phis(i_front_steady)
          phiw_N_steady = phiw(i_front_steady)
        endif

      endif
    endif

    !*** dense (H) & deposit (D) ***
    do i = i_start, mx
      hD(i) = hD(i) + dhD(i)
      z_b(i) = z_b0(i) + hD(i)
      hD_from_H(i) = hD_from_H(i) + dhD_from_H(i)
    enddo

    if(no_dense_layer==0) then
      do i = i_start, mx
        z_c(i) = z_b(i) + hH(i)
      enddo
      do i = mx, i_start, -1
        if(hH(i)>0.d0) exit
      enddo
      i_front_H = i
      x_NH = dble(i_front_H) * dx + left_boundary_point
      if(x_NH > x_NH_max) then
        x_NH_max = x_NH
      endif
      if(x_NH == x_NH_old) then
        x_NH_steady = x_NH
        i_front_H_steady = i_front_H
      endif

    elseif(no_dense_layer==1) then
      do i = i_start, mx
        z_c(i) = z_b(i)
      enddo
    endif

    do i = mx, i_start, -1
      if(hD(i)>0.d0) exit
    enddo
    i_front_D = i
    x_ND = dble(i_front_D) * dx + left_boundary_point



    !*********************************************************

    if(FC>=mx) then
      write(*,*)'The dilute layer reached the right boundary'
      exit
    endif
    if(i_front_H>=mx) then
      write(*,*)'The dense layer reached the right boundary'
      exit
    endif

    if(output_type==0) then
      if(non_dim_output==0) then
        if(time>=time_out) then
          call output(iframe,n,time,dt)
          time_out = time_out + dt_out
          iframe = iframe + 1
        endif
        if(time>=time_out_front) then
          call output_front(n,time)
          call output_mass(n,time)
          if(output_point==1) then
            call output_points(n,time,i_output_1,i_output_2,i_output_3,i_output_4)
          endif
          time_out_front = time_out_front + dt_out_front
        endif
      elseif(non_dim_output==1) then
        if(time/T_char>=time_out) then
          call output(iframe,n,time,dt)
          time_out = time_out + dt_out
          iframe = iframe + 1
        endif
        if(time/T_char>=time_out_front) then
          call output_front(n,time)
          call output_mass(n,time)
          if(output_point==1) then
            call output_points(n,time,i_output_1,i_output_2,i_output_3,i_output_4)
          endif
          time_out_front = time_out_front + dt_out_front
        endif
      else
        write(*,*)'*** ERROR *** (non_dim_output=???)'
        stop
      endif
    elseif(output_type==1) then
      call output(iframe,n,time,dt)
      iframe = iframe + 1
      call output_front(n,time)
      call output_mass(n,time)
      if(n>=10) then
        write(*,*)'************************************'
        write(*,*)'*** The calculation is finished. ***'
        write(*,*)'************************************'
        write(*,*)'   (n=',n,')'
        exit
      endif
    elseif(output_type==2) then
      call output(iframe,n,time,dt)
      iframe = iframe + 1
      call output_front(n,time)
      call output_mass(n,time)
      if(n>=100) then
        write(*,*)'************************************'
        write(*,*)'*** The calculation is finished. ***'
        write(*,*)'************************************'
        write(*,*)'   (n=',n,')'
        exit
      endif
    else
      write(*,*)'*** ERROR *** (output_type=???)'
      exit
    endif

    if(non_dim_output==0) then
      if(time > tfinal) then
        write(*,*)'************************************'
        write(*,*)'*** The calculation is finished. ***'
        write(*,*)'************************************'
        write(*,*)'   (time=',time,')'
        exit
      endif
    elseif(non_dim_output==1) then
      if(time/T_char > tfinal) then
        write(*,*)'************************************'
        write(*,*)'*** The calculation is finished. ***'
        write(*,*)'************************************'
        write(*,*)'   (time/T_char=',time/T_char,')'
        exit
      endif
    else
      write(*,*)'*** ERROR *** (non_dim_output=???)'
      stop
    endif

    if(lift_off==1 .and. no_dense_layer==1) then
      write(*,*)'*** lift_off=1 & no_dense_layer=1 ***'
      write(*,*)'   (time=',time,')'
      exit
    endif

  enddo

  close(60)
  close(70)
  close(200)
  if(output_point==1) then
    close(600)
    close(650)
    close(700)
    close(750)
    close(800)
    close(850)
    close(900)
    close(950)
  endif

  call output_runout()


  stop
end program main
