subroutine cons_to_prim(Q,h,rho,u,na,ns,nw,phia,phis,phiw,aaa,Fr,T,enthal,Cp,&
                        u_old,T_old,enthal_old,dt,i_max,step,theta,T_boil)

  use parameter
  implicit none

  integer, intent(in) :: i_max,step
  double precision, intent(in) :: dt
  double precision, dimension(0:mx+1), intent(in) :: u_old, T_old
  double precision, dimension(0:mx+1), intent(in) :: enthal_old, theta
  double precision, dimension(0:mx+1,1:6), intent(in) :: Q

  double precision, dimension(0:mx+1), intent(out) :: h, rho, u
  double precision, dimension(0:mx+1), intent(out) :: na, ns, nw, phia, phis, phiw
  double precision, dimension(0:mx+1), intent(out) :: aaa, Fr, T, enthal, Cp
  double precision, dimension(0:mx+1), intent(out) :: T_boil

  integer :: i
  double precision :: p_v


  do i = i_start, i_max
    if(Q(i,1)<=0.d0) then
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
      enthal(i) = 0.d0
      Cp(i) = 0.d0
    else
      na(i) = Q(i,2)/Q(i,1)
      ns(i) = Q(i,3)/Q(i,1)
      nw(i) = Q(i,6)/Q(i,1)
      Cp(i) = ns(i)*C_s + nw(i)*C_pw + na(i)*C_pa + (1.d0-ns(i)-nw(i)-na(i))*C_pv
      if(thermal_type==1) then
        enthal(i) = Q(i,5)/Q(i,1)
        if(step==1) then
          T(i) = T_old(i) + (enthal(i)-enthal_old(i)) / Cp(i)
        elseif(step==2) then
          T(i) = T_old(i) - dt/dx * u_old(i) * (T_old(i)-T_old(i-1)) &
                 + (enthal(i)-enthal_old(i)) / Cp(i) &
                 + dt/dx * u_old(i) * (enthal_old(i)-enthal_old(i-1)) / Cp(i)
        endif
      elseif(thermal_type==2) then
        enthal(i) = Q(i,5)/Q(i,1)
        T(i) = enthal(i)/Cp(i)
      else
        T(i) = T_0
      endif

      !*** Lift-off region ***
      if(ns(i)<=0.d0) then
        ns(i) = 0.d0
        nw(i) = 0.d0
        na(i) = 1.d0-a_v
        Cp(i) = ns(i)*C_s + nw(i)*C_pw + na(i)*C_pa + (1.d0-ns(i)-nw(i)-na(i))*C_pv
        T(i) = T_a
        rho(i) = rho_a
      endif


      p_v = p_a0 * (1.d0-ns(i)-nw(i)-na(i))*R_v
      p_v = p_v / ( (1.d0-ns(i)-nw(i)-na(i))*R_v + na(i)*R_a )
      if(1.d0-ns(i)-nw(i)-na(i) <= 0.d0) then
        p_v = 0.d0
        T_boil(i) = 0.d0
      endif


      u(i) = Q(i,4)/Q(i,1)
      rho(i) = ns(i)/rho_s + nw(i)/rho_w + ( na(i)*R_a + (1.d0-ns(i)-nw(i)-na(i))*R_v ) * T(i) / p_a0
      rho(i) = 1.d0 / rho(i)
      if(rho(i)*0.d0/=0.d0) then
        write(*,*)'*** cons_to_prim ***'
        write(*,*)'i =', i
        write(*,*)'rho(i) =', rho(i)
        write(*,*)'na(i) =', na(i)
        write(*,*)'ns(i) =', ns(i)
        write(*,*)'nw(i) =', nw(i)
        write(*,*)'Q(i,1) =', Q(i,1)
        write(*,*)'Q(i,2) =', Q(i,2)
        write(*,*)'Q(i,3) =', Q(i,3)
        write(*,*)'Q(i,6) =', Q(i,6)
      endif
      phia(i) = na(i)*rho(i) * R_a*T(i)/p_a0
      phis(i) = ns(i)*rho(i) / rho_s
      phiw(i) = nw(i)*rho(i) / rho_w
      h(i) = Q(i,1)/rho(i)
      if(rho(i)<=rho_a) then
        aaa(i) = 0.d0
      else
        aaa(i) = dsqrt( h(i) * grav * dcos(theta(i)) * (rho(i)-rho_a) / rho(i) )
      endif
      if(aaa(i)<=0.d0) then
        Fr(i) = 0.d0
      else
        Fr(i) = dabs(u(i)) / aaa(i)
      endif

    endif

  enddo


  return
end subroutine cons_to_prim
