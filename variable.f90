module variable

  implicit none

  !*** common ***
  double precision, dimension(:), allocatable :: x, z_c, theta
  

  !*** dilute ***
  integer :: FC, i_front_L, lift_off, lift_off_FC
  integer :: i_liftoff, FC_liftoff
  integer :: i_front_steady
  double precision :: x_N, dx_FC, dx_FC_new
  double precision :: x_N_max, x_N_steady, x_N_old, x_N_old_out_front
  double precision :: h_N_steady, rho_N_steady, u_N_steady
  double precision :: na_N_steady, ns_N_steady, nw_N_steady, T_N_steady
  double precision :: Fr_N_steady, enthal_N_steady
  double precision :: phia_N_steady, phis_N_steady, phiw_N_steady
  double precision :: M_0
  double precision :: time_steady
  double precision :: M_dot_top, M_dot_base, M_dot_top_H, M_dot_base_H
  double precision, dimension(:,:), allocatable :: Q, Q_star, F
  double precision, dimension(:), allocatable :: h, h_star
  double precision, dimension(:), allocatable :: rho, rho_star
  double precision, dimension(:), allocatable :: u, u_star
  double precision, dimension(:), allocatable :: na, na_star
  double precision, dimension(:), allocatable :: ns, ns_star
  double precision, dimension(:), allocatable :: nw, nw_star
  double precision, dimension(:), allocatable :: phia, phia_star
  double precision, dimension(:), allocatable :: phis, phis_star
  double precision, dimension(:), allocatable :: phiw, phiw_star
  double precision, dimension(:), allocatable :: aaa, aaa_star
  double precision, dimension(:), allocatable :: Fr, Fr_star
  double precision, dimension(:), allocatable :: T, T_star
  double precision, dimension(:), allocatable :: enthal, enthal_star
  double precision, dimension(:), allocatable :: Cp, Cp_star, Cp_old
  double precision, dimension(:), allocatable :: T_boil, T_boil_star
  double precision, dimension(:), allocatable :: c

  !*** dense ***
  integer :: i_front_H, i_front_H_steady
  double precision :: x_NH, x_NH_max, x_NH_steady, x_NH_old, x_NH_old_out_front
  double precision :: time_H_steady
  double precision, dimension(:,:), allocatable :: QH, QH_star, FH
  double precision, dimension(:), allocatable :: hH, hH_star
  double precision, dimension(:), allocatable :: uH, uH_star
  double precision, dimension(:), allocatable :: aaaH, aaaH_star
  double precision, dimension(:), allocatable :: FrH, FrH_star
  double precision, dimension(:,:), allocatable :: interact

  !*** deposit ***
  integer :: i_front_D
  double precision :: x_ND
  double precision, dimension(:), allocatable :: z_b, z_b0
  double precision, dimension(:), allocatable :: hD, dhD
  double precision, dimension(:), allocatable :: hD_from_H, dhD_from_H

end module variable
