subroutine not_kinematic_condition(WW,rho,dt,A1,A4,theta)

  use parameter
  implicit none

  double precision, intent(in) :: rho,dt,A1,A4,theta
  double precision, dimension(1:3), intent(inout) :: WW !=(H,U,X)

  integer, parameter :: ntrial=1000
  integer :: n, i, info
  double precision, parameter :: TOL_f=1.d-10, TOL_w=1.d-10
  double precision :: err_f, err_w
  double precision :: H, U, X
  double precision, dimension(1:3,1:3) :: AA
  double precision, dimension(1:3) :: b, piv

  do n = 1, ntrial
    H = WW(1)
    U = WW(2)
    X = WW(3)
!    write(*,*)'*** Newton-Raphson ***  n=',n
!    write(*,*)'H=',H
!    write(*,*)'U=',U
!    write(*,*)'X=',X

    b(1) = - ( rho*H*X - A1 )
    b(2) = - ( U*U - (Fr_N0*Fr_N0*grav*dcos(theta)/rho_a)*(rho-rho_a)*H )
    b(3) = - ( rho*H*U*X - A4 + 0.5d0*dt*grav*dcos(theta)*(rho-rho_a)*H*H )
!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_f = 0.d0
    do i = 1, 3
      err_f = err_f + dabs(b(i))
    enddo
    if(err_f<=TOL_f) return

    !*** the function of mass conservation equation ***
    AA(1,1) = rho * X ! df/dH
    AA(1,2) = 0.d0 ! df/dU
    AA(1,3) = rho * H ! df/dX

    !*** the function of front condition ***
    AA(2,1) = - (Fr_N0*Fr_N0*grav*dcos(theta)/rho_a)*(rho-rho_a) ! df/dH
    AA(2,2) = 2.d0*U ! df/dU
    AA(2,3) = 0.d0 ! df/dX

    !*** the function of momentun conservation equation ***
    AA(3,1) = rho*U*X + dt*grav*dcos(theta)*(rho-rho_a)*H ! df/dH
    AA(3,2) = rho * H * X ! df/dU
    AA(3,3) = rho * H * U ! df/dX

    !*** using LAPACK & BLAS ***
    call dgesv(3,1,AA,3,piv,b,3,info)
!    write(*,*)'info=',info

!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_w = 0.d0
    do i = 1, 3
      err_w = err_w + dabs(b(i))
      WW(i) = WW(i) + b(i)
    enddo
    if(err_w<=TOL_w) return

!    write(*,*)'  err_w=',err_w
!    write(*,*)'  err_f=',err_f

    if(n>=ntrial) then
      write(*,*)'*** ERROR *** (Newton-Raphson)'
      write(*,*)'  n=',n,'  / ntrial=',ntrial
      write(*,*)'  err_w=',err_w
      write(*,*)'  err_f=',err_f
      write(*,*)'  H=',WW(1)
      write(*,*)'  U=',WW(2)
      write(*,*)'  X=',WW(3)
      write(*,*)'  rho=',rho
      write(*,*)'  A1=',A1
      write(*,*)'  A4=',A4
      stop
    endif


  enddo


  return
end subroutine not_kinematic_condition



subroutine not_momentum_equation(WW,rho,dt,A1,dx_FC,theta)

  use parameter
  implicit none

  double precision, intent(in) :: rho,dt,A1,dx_FC,theta
  double precision, dimension(1:3), intent(inout) :: WW !=(H,U,X)

  integer, parameter :: ntrial=1000
  integer :: n, i, info
  double precision, parameter :: TOL_f=1.d-10, TOL_w=1.d-10
  double precision :: err_f, err_w
  double precision :: H, U, X
  double precision, dimension(1:3,1:3) :: AA
  double precision, dimension(1:3) :: b, piv

  do n = 1, ntrial
    H = WW(1)
    U = WW(2)
    X = WW(3)
!    write(*,*)'*** Newton-Raphson ***  n=',n
!    write(*,*)'H=',H
!    write(*,*)'U=',U
!    write(*,*)'X=',X

    b(1) = - (rho*H*X - A1)
    b(2) = - (U*U - (Fr_N0*Fr_N0*grav*dcos(theta)/rho_a) * (rho-rho_a)*H)
    b(3) = - (X - dt*U - dx_FC)
!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_f = 0.d0
    do i = 1, 3
      err_f = err_f + dabs(b(i))
    enddo
    if(err_f<=TOL_f) return

    !*** the function of mass conservation equation ***
    AA(1,1) = rho * X ! df/dH
    AA(1,2) = 0.d0 ! df/dU
    AA(1,3) = rho * H ! df/dX

    !*** the function of front condition ***
    AA(2,1) = - (Fr_N0*Fr_N0*grav*dcos(theta)/rho_a)*(rho-rho_a) ! df/dH
    AA(2,2) = 2.d0*U ! df/dU
    AA(2,3) = 0.d0 ! df/dX

    !*** the function of kinematic condition ***
    AA(3,1) = 0.d0
    AA(3,2) = - dt
    AA(3,3) = 1.d0

    !*** using LAPACK & BLAS ***
    call dgesv(3,1,AA,3,piv,b,3,info)
!    write(*,*)'info=',info

!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_w = 0.d0
    do i = 1, 3
      err_w = err_w + dabs(b(i))
      WW(i) = WW(i) + b(i)
    enddo
    if(err_w<=TOL_w) return

!    write(*,*)'  err_w=',err_w
!    write(*,*)'  err_f=',err_f

    if(n>=ntrial) then
      write(*,*)'*** ERROR *** (Newton-Raphson)'
      write(*,*)'  n=',n,'  / ntrial=',ntrial
      write(*,*)'  err_w=',err_w
      write(*,*)'  err_f=',err_f
      write(*,*)'  H=',WW(1)
      write(*,*)'  U=',WW(2)
      write(*,*)'  X=',WW(3)
      write(*,*)'  rho=',rho
      write(*,*)'  A1=',A1
      write(*,*)'  dx_FC=',dx_FC
      stop
    endif


  enddo


  return
end subroutine not_momentum_equation



subroutine not_front_condition(WW,rho,dt,A1,A4,dx_FC,theta)

  use parameter
  implicit none

  double precision, intent(in) :: rho,dt,A1,A4,dx_FC,theta
  double precision, dimension(1:3), intent(inout) :: WW !=(H,U,X)

  integer, parameter :: ntrial=1000
  integer :: n, i, info
  double precision, parameter :: TOL_f=1.d-10, TOL_w=1.d-10
  double precision :: err_f, err_w
  double precision :: H, U, X
  double precision, dimension(1:3,1:3) :: AA
  double precision, dimension(1:3) :: b, piv

  do n = 1, ntrial
    H = WW(1)
    U = WW(2)
    X = WW(3)
!    write(*,*)'*** Newton-Raphson ***  n=',n
!    write(*,*)'H=',H
!    write(*,*)'U=',U
!    write(*,*)'X=',X

    b(1) = - (rho*H*X - A1)
    b(2) = - (rho*H*U*X - A4 + 0.5d0*dt*grav*dcos(theta)*(rho-rho_a)*H*H)
    b(3) = - (X - dt*U - dx_FC)
!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_f = 0.d0
    do i = 1, 3
      err_f = err_f + dabs(b(i))
    enddo
    if(err_f<=TOL_f) return

    !*** the function of mass conservation equation ***
    AA(1,1) = rho * X ! df/dH
    AA(1,2) = 0.d0 ! df/dU
    AA(1,3) = rho * H ! df/dX

    !*** the function of momentun conservation equation ***
    AA(2,1) = rho*U*X + dt*grav*dcos(theta)*(rho-rho_a)*H ! df/dH
    AA(2,2) = rho * H * X ! df/dU
    AA(2,3) = rho * H * U ! df/dX

    !*** the function of kinematic condition ***
    AA(3,1) = 0.d0
    AA(3,2) = - dt
    AA(3,3) = 1.d0

    !*** using LAPACK & BLAS ***
    call dgesv(3,1,AA,3,piv,b,3,info)
!    write(*,*)'info=',info

!    write(*,*)'b(1)=',b(1)
!    write(*,*)'b(2)=',b(2)
!    write(*,*)'b(3)=',b(3)

    err_w = 0.d0
    do i = 1, 3
      err_w = err_w + dabs(b(i))
      WW(i) = WW(i) + b(i)
    enddo
    if(err_w<=TOL_w) return

!    write(*,*)'  err_w=',err_w
!    write(*,*)'  err_f=',err_f

    if(n>=ntrial) then
      write(*,*)'*** ERROR *** (Newton-Raphson)'
      write(*,*)'  n=',n,'  / ntrial=',ntrial
      write(*,*)'  err_w=',err_w
      write(*,*)'  err_f=',err_f
      write(*,*)'  H=',WW(1)
      write(*,*)'  U=',WW(2)
      write(*,*)'  X=',WW(3)
      write(*,*)'  rho=',rho
      write(*,*)'  A1=',A1
      write(*,*)'  A4=',A4
      write(*,*)'  dx_FC=',dx_FC
      stop
    endif


  enddo


  return
end subroutine not_front_condition

