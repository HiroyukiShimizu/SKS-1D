subroutine output_mass(n,time)

  use parameter
  use variable
  implicit none

  integer, intent(in) :: n
  double precision, intent(in) :: time
  integer :: i, i_max
  double precision :: M_tot, M_HD, M_D, M_D_from_H, M_D_from_H_overshoot
  double precision :: hD_ave, aspect_ratio
  double precision :: Ws, Ent, D

  if(time==0.d0) then
    write(200,*)'# time, M_tot, M_HD, M_D, M_HD/M_tot, M_D/M_tot, hD_ave, aspect_ratio,'
    write(200,*)'# M_D_from_H, M_D_from_H_os, M_D_from_H/M_tot, M_D_from_H_os/M_tot'
    write(200,*)'# M_dot_top, M_dot_base, M_dot_top_H, M_dot_base_H'
    write(200,*)'# M_dot_top/M_dot_0, M_dot_base/M_dot_0, M_dot_top_H/M_dot_0, M_dot_base_H/M_dot_0'
  endif

  i_max = max(i_front_H, i_front_D)
  M_HD = 0.d0
  M_D = 0.d0
  M_D_from_H = 0.d0
  M_D_from_H_overshoot = 0.d0
  hD_ave = 0.d0
  aspect_ratio = 0.d0
  M_dot_top = 0.d0
  M_dot_base = 0.d0
  M_dot_top_H = 0.d0
  M_dot_base_H = 0.d0

  if(configuration_type==1) then
    if(time==0.d0) then
      M_0 = 0.d0
      do i = i_start, i_front_L
        M_0 = M_0 + rho(i) * h(i)
      enddo
      M_0 = M_0 * y_0 * dx
    endif
    M_tot = M_0 + M_dot_0 * time
    if(i_start > i_max) then
      M_HD = 0.d0
      M_D = 0.d0
      M_D_from_H = 0.d0
      M_D_from_H_overshoot = 0.d0
    else
      do i = i_start, i_max
        M_HD = M_HD + phi_sH*hH(i) + phi_sD*hD(i)
      enddo
      M_HD = M_HD * y_0 * rho_s * dx
      do i = i_start, i_front_D
        M_D = M_D + hD(i)
        M_D_from_H = M_D_from_H + hD_from_H(i)
        if(x(i)>x_NH) then
          M_D_from_H_overshoot = M_D_from_H_overshoot + hD_from_H(i)
        endif
      enddo
      M_D = M_D * y_0 * phi_sD*rho_s * dx
      M_D_from_H = M_D_from_H * y_0 * phi_sD*rho_s * dx
      M_D_from_H_overshoot = M_D_from_H_overshoot * y_0 * phi_sD*rho_s * dx
      hD_ave = M_D / (y_0*phi_sD*rho_s*(x(i_front_D)-x_0))
      aspect_ratio = hD_ave / (x(i_front_D)-x_0)
    endif

    if(i_start < i_front_H) then
      do i = i_start, i_front_H
        if(hH(i)>h_Hchar*eps) then
          M_dot_base_H = M_dot_base_H + dx
          if(rho(i)>rho_a) then
            Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
            M_dot_top_H = M_dot_top_H + ns(i)*rho(i)*Ws*dcos(theta(i))*dx
          endif
        endif
      enddo
      D = D_func(rho_H,T_0,phi_sH*rho_s/rho_H,0.d0,0.d0)
      M_dot_base_H = M_dot_base_H * y_0*phi_sD*rho_s*D*dcos(theta(i))
      M_dot_top_H = M_dot_top_H * y_0
    endif

    do i = i_start, i_front_L
      if(rho(i)>rho_a) then
        Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
        M_dot_base = M_dot_base + ns(i)*rho(i)*Ws*dcos(theta(i))*dx
        Ent = Ent_func(Fr(i))
        M_dot_top = M_dot_top + dabs(u(i))*Ent*dx
      endif
    enddo
    M_dot_base = M_dot_base * y_0
    M_dot_top = M_dot_top * y_0*rho_a

  elseif(configuration_type==2) then
    if(time==0.d0) then
      M_0 = 0.d0
      do i = i_start, i_front_L
        M_0 = M_0 + rho(i) * h(i) * ( x(i)*x(i) - x(i-1)*x(i-1) )
      enddo
      M_0 = M_0 * pi
    endif
    M_tot = M_0 + M_dot_0 * time
    if(i_start > i_max) then
      M_HD = 0.d0
      M_D = 0.d0
      M_D_from_H = 0.d0
      M_D_from_H_overshoot = 0.d0
    else
      do i = i_start, i_max
        M_HD = M_HD + ( phi_sH*hH(i) + phi_sD*hD(i) ) * ( x(i)*x(i) - x(i-1)*x(i-1) )
      enddo
      M_HD = M_HD * rho_s * pi
      do i = i_start, i_front_D
        M_D = M_D + hD(i) * ( x(i)*x(i) - x(i-1)*x(i-1) )
        M_D_from_H = M_D_from_H + hD_from_H(i) * ( x(i)*x(i) - x(i-1)*x(i-1) )
        if(x(i)>x_NH) then
          M_D_from_H_overshoot = M_D_from_H_overshoot + hD_from_H(i) * ( x(i)*x(i) - x(i-1)*x(i-1) )
        endif
      enddo
      M_D = M_D * phi_sD*rho_s * pi
      M_D_from_H = M_D_from_H * phi_sD*rho_s * pi
      M_D_from_H_overshoot = M_D_from_H_overshoot * phi_sD*rho_s * pi
      hD_ave = M_D / (pi*phi_sD*rho_s*(x(i_front_D)*x(i_front_D)-x_0*x_0))
      aspect_ratio = hD_ave / (2.d0*x(i_front_D))
    endif

    if(i_start < i_front_H) then
      do i = i_start, i_front_H
        if(hH(i)>h_Hchar*eps) then
          M_dot_base_H = M_dot_base_H + x(i)*dx
          if(rho(i)>rho_a.and.i<=i_front_L) then
            Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
            M_dot_top_H = M_dot_top_H + ns(i)*rho(i)*Ws*dcos(theta(i))*x(i)*dx
          endif
        endif
      enddo
      D = D_func(rho_H,T_0,phi_sH*rho_s/rho_H,0.d0,0.d0)
      M_dot_base_H = M_dot_base_H * 2.d0*pi*phi_sD*rho_s*D*dcos(theta(i))
      M_dot_top_H = M_dot_top_H * 2.d0*pi
    endif

    do i = i_start, i_front_L
      if(rho(i)>rho_a) then
        Ws = Ws_func(rho(i),T(i),ns(i),na(i),nw(i))
        M_dot_base = M_dot_base + ns(i)*rho(i)*Ws*dcos(theta(i))*x(i)*dx
        Ent = Ent_func(Fr(i))
        M_dot_top = M_dot_top + dabs(u(i))*Ent*x(i)*dx
      endif
    enddo
    M_dot_base = M_dot_base * 2.d0*pi
    M_dot_top = M_dot_top * 2.d0*pi*rho_a

  else
    write(*,*)'*** ERROR *** (in output_mass.f90)'
    write(*,*)'  configuration_type =', configuration_type
    stop
  endif

  if(non_dim_output==0) then
    write(200,'(20E26.16)') time, M_tot, M_HD, M_D, M_HD/M_tot, M_D/M_tot, hD_ave, aspect_ratio, &
                            M_D_from_H, M_D_from_H_overshoot, M_D_from_H/M_tot, M_D_from_H_overshoot/M_tot, &
                            M_dot_top, M_dot_base, M_dot_top_H, M_dot_base_H, &
                            M_dot_top/M_dot_0, M_dot_base/M_dot_0, M_dot_top_H/M_dot_0, M_dot_base_H/M_dot_0
                            
  elseif(non_dim_output==1) then
    write(200,'(20E26.16)') time/T_char, M_tot, M_HD, M_D, M_HD/M_tot, M_D/M_tot, hD_ave/H_char, aspect_ratio, &
                            M_D_from_H, M_D_from_H_overshoot, M_D_from_H/M_tot, M_D_from_H_overshoot/M_tot, &
                            M_dot_top, M_dot_base, M_dot_top_H, M_dot_base_H, &
                            M_dot_top/M_dot_0, M_dot_base/M_dot_0, M_dot_top_H/M_dot_0, M_dot_base_H/M_dot_0
  endif


  return
end subroutine output_mass
