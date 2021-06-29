subroutine hllc_flux_H(FH,QH,theta)

  use parameter
  implicit none

  double precision, dimension(0:mx+1,1:2), intent(in) :: QH
  double precision, dimension(0:mx+1), intent(in) :: theta
  double precision, dimension(0:mx+1,1:2), intent(out) :: FH

  integer :: i
  double precision :: h_L, u_L, aaa_L
  double precision :: h_R, u_R, aaa_R
  double precision :: q1_L, q2_L
  double precision :: q1_R, q2_R
  double precision :: f1_L, f2_L
  double precision :: f1_R, f2_R
  double precision :: S_L, S_R
  double precision :: h_star0, q_L, q_R
  double precision :: theta_L, theta_R
  !*** HLLC ***
  double precision :: h_star, u_star, S_star
  double precision :: f1_sL, f2_sL
  double precision :: f1_sR, f2_sR


  do i = i_start, mx+1
    !*** conserved variables ***
    q1_L = QH(i-1,1)
    q2_L = QH(i-1,2)
    theta_L = theta(i-1)

    q1_R = QH(i,1)
    q2_R = QH(i,2)
    theta_R = theta(i)

    !*** primitive variables ***
    h_L = q1_L
    u_L = q2_L/h_L
    aaa_L = dsqrt(g_H*h_L*dcos(theta_L))

    h_R = q1_R
    u_R = q2_R/h_R
    aaa_R = dsqrt(g_H*h_R*dcos(theta_R))


    !*** flux ***
    f1_L = u_L * h_L
    f2_L = u_L*u_L*h_L + 0.5d0*g_H*h_L*h_L*dcos(theta_L)

    f1_R = u_R * h_R
    f2_R = u_R*u_R*h_R + 0.5d0*g_H*h_R*h_R*dcos(theta_R)

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
    if(h_star>h_R) then
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
    h_star = h_R * (S_R-u_R)/(S_R-u_star)
    S_star = u_star

    !*** HLLC flux ***
    f1_sL = u_star * h_star
    f2_sL = u_star*u_star*h_star + 0.5d0*g_H*h_star*h_star*dcos(theta_L)

    f1_sR = u_star * h_star
    f2_sR = u_star*u_star*h_star + 0.5d0*g_H*h_star*h_star*dcos(theta_R)


    if(S_L>=0.d0) then
      FH(i,1) = f1_L
      FH(i,2) = f2_L
    elseif(S_R<=0.d0) then
      FH(i,1) = f1_R
      FH(i,2) = f2_R
    elseif(S_L<0.d0.and.0.d0<=S_star) then
      FH(i,1) = f1_sL
      FH(i,2) = f2_sL
    elseif(S_star<0.d0.and.0.d0<S_R) then
      FH(i,1) = f1_sR
      FH(i,2) = f2_sR
    else
      write(*,*)'*** ERROR *** (HLLC in dense layer)'
      write(*,*)'  S_L =', S_L
      write(*,*)'  S_R =', S_R
      write(*,*)'  S_star =', S_star
      stop
    endif

  enddo

  return
end subroutine hllc_flux_H
