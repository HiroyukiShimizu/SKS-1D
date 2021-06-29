subroutine hll_flux_L(i_max,F,h,rho,u,na,ns,nw,aaa,enthal,theta)

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
  double precision :: h_star, q_L, q_R
  double precision :: theta_L, theta_R


!  do i = 1, i_max
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
    if(aaa_L*0.d0/=0.d0) aaa_L=0.d0

    h_R = h(i)
    rho_R = rho(i)
    u_R = u(i)
    na_R = na(i)
    ns_R = ns(i)
    nw_R = nw(i)
    aaa_R = aaa(i)
    enthal_R = enthal(i)
    theta_R = theta(i)
    if(aaa_R*0.d0/=0.d0) aaa_R=0.d0

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
    h_star = 0.5d0*(h_L+h_R)-0.25d0*(u_R-u_L)*(h_L+h_R)/(aaa_L+aaa_R)
!    if(h_star*0.d0/=0.d0) h_star=0.d0
    if(h_star>h_L) then
      q_L = 0.5d0*(h_star+h_L)*h_star/h_L/h_L
      q_L = dsqrt(q_L)
      if(q_L*0.d0/=0.d0) q_L = 1.d0
    else
      q_L = 1.d0
    endif
    if(h_star>h_R) then
      q_R = 0.5d0*(h_star+h_R)*h_star/h_R/h_R
      q_R = dsqrt(q_R)
    else
      q_R = 1.d0
    endif
    S_L = u_L-aaa_L*q_L
    S_R = u_R+aaa_R*q_R

    if(S_R==0.d0.and.S_L==0.d0) then
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
    else
!      if(S_R-S_L==0.d0) then
!        F(i,1) = 0.d0
!        F(i,2) = 0.d0
!        F(i,3) = 0.d0
!        F(i,4) = 0.d0
!        F(i,5) = 0.d0
!      else
        F(i,1) = ( S_R*f1_L - S_L*f1_R + S_L*S_R*(q1_R-q1_L) ) / (S_R-S_L)
        F(i,2) = ( S_R*f2_L - S_L*f2_R + S_L*S_R*(q2_R-q2_L) ) / (S_R-S_L)
        F(i,3) = ( S_R*f3_L - S_L*f3_R + S_L*S_R*(q3_R-q3_L) ) / (S_R-S_L)
        F(i,4) = ( S_R*f4_L - S_L*f4_R + S_L*S_R*(q4_R-q4_L) ) / (S_R-S_L)
        F(i,5) = ( S_R*f5_L - S_L*f5_R + S_L*S_R*(q5_R-q5_L) ) / (S_R-S_L)
        F(i,6) = ( S_R*f6_L - S_L*f6_R + S_L*S_R*(q6_R-q6_L) ) / (S_R-S_L)
!      endif

      if(F(i,1)*0.d0/=0.d0) then
        write(*,*)'*** ERROR *** (in hll_flux_L.f90)'
        write(*,*)'i =', i
        write(*,*)'F(i,1) =', F(i,1)
        write(*,*)'rho_L =', rho_L
        write(*,*)'rho_R =', rho_R
        write(*,*)'S_L =', S_L
        write(*,*)'S_R =', S_R
        write(*,*)'u_L =', u_L
        write(*,*)'u_R =', u_R
        write(*,*)'h_L =', h_L
        write(*,*)'h_R =', h_R
        write(*,*)'aaa_L =', aaa_L
        write(*,*)'aaa_R =', aaa_R
        write(*,*)'q_L =', q_L
        write(*,*)'q_R =', q_R
        write(*,*)'h_star =', h_star
        write(*,*)''
      endif
    endif

  enddo

  return
end subroutine hll_flux_L
