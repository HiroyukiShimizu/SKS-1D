subroutine cfl_condition(dt)

  use parameter
  use variable
  implicit none

  double precision, intent(out) :: dt

  integer :: i
  double precision :: lmd_max, lmd_1, lmd_2, lmd_3

  lmd_max = 0.d0

  if(lift_off==0) then
    do i = 1, FC
      lmd_1 = dabs( u(i) - aaa(i) )
      lmd_3 = dabs( u(i) + aaa(i) )
      if(lmd_max<lmd_1) lmd_max = lmd_1
      if(lmd_max<lmd_3) lmd_max = lmd_3
    enddo
  endif

  if(no_dense_layer==0) then
    do i = 1, i_front_H
      lmd_1 = dabs( uH(i) - dsqrt(QH(i,1)*g_H*dcos(theta(i))) )
      lmd_2 = dabs( uH(i) + dsqrt(QH(i,1)*g_H*dcos(theta(i))) )
      if(lmd_max<lmd_1) lmd_max = lmd_1
      if(lmd_max<lmd_2) lmd_max = lmd_2
    enddo
  endif

  dt = cfl * dx / lmd_max
  if(non_dim_output==0) then
    if(dt>dt_out) then
      dt = 0.1d0*dt_out
    endif
  else
    if(dt>dt_out*T_char) then
      dt = 0.1d0*dt_out*T_char
    endif
  endif


  return
end subroutine cfl_condition
