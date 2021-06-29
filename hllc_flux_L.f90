subroutine hllc_flux_L(i_max,F,h,rho,u,na,ns,nw,aaa,enthal,theta)

  use parameter
  implicit none

  integer, intent(in) :: i_max
  double precision, dimension(0:mx+1), intent(in) :: h,rho,u,na,ns,nw,aaa,enthal
  double precision, dimension(0:mx+1), intent(in) :: theta
  double precision, dimension(0:mx+1,1:6), intent(out) :: F

  integer :: i
  double precision :: h_L, rho_L, u_L, na_L, ns_L, nw_L, aaa_L, enthal_L
  double precision :: h_R, rho_R, u_R, na_R, ns_R, nw_R, aaa_R, enthal_R
  double precision :: q1_L, q2_L, q3_L, q4_L, q5_L, q6_L
  double precision :: q1_R, q2_R, q3_R, q4_R, q5_R, q6_R
  double precision :: f1_L, f2_L, f3_L, f4_L, f5_L, f6_L
  double precision :: f1_R, f2_R, f3_R, f4_R, f5_R, f6_R
  double precision :: S_L, S_R
  double precision :: h_star0, q_L, q_R
  double precision :: theta_L, theta_R
  !*** HLLC ***
  double precision :: h_star, u_star, S_star
  double precision :: f1_sL, f2_sL, f3_sL, f4_sL, f5_sL, f6_sL
  double precision :: f1_sR, f2_sR, f3_sR, f4_sR, f5_sR, f6_sR


  do i = i_start, i_max
    !*** copy ***
    h_L = h(i-1)
    rho_L = rho(i-1)
    u_L = u(i-1)
    na_L = na(i-1)
    ns_L = ns(i-1)
    nw_L = nw(i-1)
    aaa_L = aaa(i-1)
    enthal_L = enthal(i-1)
    theta_L = theta(i-1)

    h_R = h(i)
    rho_R = rho(i)
    u_R = u(i)
    na_R = na(i)
    ns_R = ns(i)
    nw_R = nw(i)
    aaa_R = aaa(i)
    enthal_R = enthal(i)
    theta_R = theta(i)

    !*** conserved variables ***
    q1_L = rho_L * h_L
    q2_L = na_L * rho_L * h_L
    q3_L = ns_L * rho_L * h_L
    q4_L = rho_L * u_L * h_L
    q5_L = rho_L * enthal_L * h_L
    q6_L = nw_L * rho_L * h_L

    q1_R = rho_R * h_R
    q2_R = na_R * rho_R * h_R
    q3_R = ns_R * rho_R * h_R
    q4_R = rho_R * u_R * h_R
    q5_R = rho_R * enthal_R * h_R
    q6_R = nw_R * rho_R * h_R

    !*** flux ***
    f1_L = rho_L * u_L * h_L
    f2_L = na_L * rho_L * u_L * h_L
    f3_L = ns_L * rho_L * u_L * h_L
    f4_L = rho_L*u_L*u_L*h_L + 0.5d0*(rho_L-rho_a)*grav*h_L*h_L*dcos(theta_L)
    f5_L = rho_L * u_L * enthal_L * h_L
    f6_L = nw_L * rho_L * u_L * h_L

    f1_R = rho_R * u_R * h_R
    f2_R = na_R * rho_R * u_R * h_R
    f3_R = ns_R * rho_R * u_R * h_R
    f4_R = rho_R*u_R*u_R*h_R + 0.5d0*(rho_R-rho_a)*grav*h_R*h_R*dcos(theta_R)
    f5_R = rho_R * u_R * enthal_R * h_R
    f6_R = nw_R * rho_R * u_R * h_R

    !*** chracteric velocity ***
!    S_L = u_L-aaa_L
!    S_R = u_R+aaa_R
!    S_L = dmin1(u_L-aaa_L, u_R-aaa_R)
!    S_R = dmax1(u_L+aaa_L, u_R+aaa_R)
    h_star0 = 0.5d0*(h_L+h_R)-0.25d0*(u_R-u_L)*(h_L+h_R)/(aaa_L+aaa_R)
    if(h_star>h_L) then
      q_L = 0.5d0*(h_star0+h_L)*h_star0/h_L/h_L
      q_L = dsqrt(q_L)
    else
      q_L = 1.d0
    endif
    if(h_star0>h_R) then
      q_R = 0.5d0*(h_star0+h_R)*h_star0/h_R/h_R
      q_R = dsqrt(q_R)
    else
      q_R = 1.d0
    endif
    S_L = u_L-aaa_L*q_L
    S_R = u_R+aaa_R*q_R

    !*** HLLC star region ***
    u_star = S_L*h_R*(u_R-S_R) - S_R*h_L*(u_L-S_L)
    u_star = u_star / (h_R*(u_R-S_R) - h_L*(u_L-S_L))
!    h_star = h_R * (S_R-u_R)/(S_R-u_star)
    h_star = h_L * (S_L-u_L)/(S_L-u_star)
    S_star = u_star

    !*** HLLC flux ***
    f1_sL = rho_L * u_star * h_star
    f2_sL = na_L * rho_L * u_star * h_star
    f3_sL = ns_L * rho_L * u_star * h_star
    f4_sL = rho_L * u_star * u_star * h_star
    f4_sL = f4_sL + 0.5d0*(rho_L-rho_a)*grav*h_star*h_star*dcos(theta_L)
    f5_sL = rho_L * u_star * enthal_L * h_star
    f6_sL = nw_L * rho_L * u_star * h_star

    f1_sR = rho_R * u_star * h_star
    f2_sR = na_R * rho_R * u_star * h_star
    f3_sR = ns_R * rho_R * u_star * h_star
    f4_sR = rho_R * u_star * u_star * h_star
    f4_sR = f4_sR + 0.5d0*(rho_R-rho_a)*grav*h_star*h_star*dcos(theta_R)
    f5_sR = rho_R * u_star * enthal_R * h_star
    f6_sR = nw_R * rho_R * u_star * h_star


    if(S_L==0.d0.and.S_R==0.d0) then
      F(i,1) = 0.d0
      F(i,2) = 0.d0
      F(i,3) = 0.d0
      F(i,4) = 0.d0
      F(i,5) = 0.d0
      F(i,6) = 0.d0
    elseif(S_L>=0.d0) then
      F(i,1) = f1_L
      F(i,2) = f2_L
      F(i,3) = f3_L
      F(i,4) = f4_L
      F(i,5) = f5_L
      F(i,6) = f6_L
    elseif(S_R<=0.d0) then
      F(i,1) = f1_R
      F(i,2) = f2_R
      F(i,3) = f3_R
      F(i,4) = f4_R
      F(i,5) = f5_R
      F(i,6) = f6_R
    elseif(S_L<0.d0.and.0.d0<=S_star) then
      F(i,1) = f1_sL
      F(i,2) = f2_sL
      F(i,3) = f3_sL
      F(i,4) = f4_sL
      F(i,5) = f5_sL
      F(i,6) = f6_sL
    elseif(S_star<0.d0.and.0.d0<S_R) then
      F(i,1) = f1_sR
      F(i,2) = f2_sR
      F(i,3) = f3_sR
      F(i,4) = f4_sR
      F(i,5) = f5_sR
      F(i,6) = f6_sR
    else
      write(*,*)'*** ERROR *** (HLLC in dilute layer)'
      write(*,*)'  S_L =', S_L
      write(*,*)'  S_R =', S_R
      write(*,*)'  S_star =', S_star
      stop
    endif

    if(F(i,1)*0.d0/=0.d0) then
      write(*,*)'i =', i
      write(*,*)'F(i,1) =', F(i,1)
      write(*,*)'rho_L =', rho_L
      write(*,*)'rho_R =', rho_R
    endif

  enddo

  return
end subroutine hllc_flux_L
