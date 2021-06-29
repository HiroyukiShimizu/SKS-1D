#! /bin/bash

horizontal_axis="x"
#horizontal_axis="r"

gnuplot << EOF
#set terminal postscript eps color enhanced "Arial" 30
#set output "h_points.eps"
#set xrange [0:10]
#set yrange [0:4]
#set ytics 0,1,4
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Thickness,  {/Arial-Italic h} [m]}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:2 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:2 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:2 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:2 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"


chai_0 = 3.6;
chai_1 = 0;
z_0 = 0.3215;
z_1 = -0.01421
func(x) = (1+chai_0+chai_1*x)*(z_0+z_1*x) 

set terminal postscript eps color enhanced "Arial" 30
set output "h_Profile-2.eps"
set xrange [0:6]
set yrange [0:4]
set ytics 0,1,4
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Thickness,  {/Arial-Italic h} [m]}'
unset title
set key
plot \
func(x) w l lt 1 lc rgb "dark-grey" lw 20 ti "Experimental", \
1.35378 w l lt 0 lc rgb "black" lw 20 ti "Exp. (time-ave.)", \
'point_1_L.dat' using 1:2 w l lt 1 lc rgb "red" lw 6 ti "Numerical"

chai_0 = 3.6;
chai_1 = 0;
z_0 = 0.2721;
z_1 = 0.003674
func(x) = (1+chai_0+chai_1*x)*(z_0+z_1*x) 

set terminal postscript eps color enhanced "Arial" 30
set output "h_Profile-3.eps"
set xrange [0:7]
set yrange [0:4]
set ytics 0,1,4
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Thickness,  {/Arial-Italic h} [m]}'
unset title
set key
plot \
func(x) w l lt 1 lc rgb "dark-grey" lw 20 ti "Experimental", \
1.2857 w l lt 0 lc rgb "black" lw 20 ti "Exp. (time-ave.)", \
'point_2_L.dat' using 1:2 w l lt 1 lc rgb "red" lw 6 ti "Numerical"

#set terminal postscript eps color enhanced "Arial" 30
#set output "rho_points.eps"
#set xrange [0:10]
#set yrange [0:3.5]
#set ytics 0,0.5,3.5
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Density,  {/Symbol-Oblique r} [kg/m^3]}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:3 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:3 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:3 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:3 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "u_points.eps"
#set xrange [0:10]
#set yrange [0:8]
#set ytics 0,1,8
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Velocity,  {/Arial-Italic u} [m/s]}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:4 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:4 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:4 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:4 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "n_a_points.eps"
#set xrange [0:10]
#set yrange [0:1]
#set ytics 0,0.1,1
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Air mass fraction,  {/Arial-Italic n}_{a}}'
#unset title
#set key top left
#plot \
#'point_1_L.dat' using 1:5 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:5 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:5 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:5 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "n_s_points.eps"
#set xrange [0:10]
#set yrange [0:1]
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Solid mass fraction,  {/Arial-Italic n}_{s}}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:6 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:6 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:6 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:6 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "T_points.eps"
#set xrange [0:10]
#set yrange [260:300]
#set ytics 260,10,300
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Temperature,  {/Arial-Italic T} [K]}'
#unset title
#set key bottom right
#plot \
#284 w l lt 0 lc rgb "black" lw 6 ti "{/Arial-Italic T}_{a}", \
#'point_1_L.dat' using 1:7 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:7 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:7 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:7 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "rho_rho_a_points.eps"
#set xrange [0:10]
#set yrange [1:2]
#set ytics 0,0.1,2
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Density ratio,  {/Symbol-Oblique r}/{/Symbol-Oblique r}_a}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:10 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:10 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:10 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:10 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "z_c_points.eps"
#set xrange [0:10]
#set yrange [0:0.01]
#set ytics 0,0.002,0.01
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Height,  {/Arial-Italic z}_c [m]}'
#unset title
#set key
#plot \
#'point_1_H.dat' using 1:6 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_H.dat' using 1:6 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_H.dat' using 1:6 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_H.dat' using 1:6 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "z_b_points.eps"
#set xrange [0:10]
#set yrange [0:0.01]
#set ytics 0,0.002,0.01
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Height,  {/Arial-Italic z}_b [m]}'
#unset title
#set key
#plot \
#'point_1_H.dat' using 1:8 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_H.dat' using 1:8 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_H.dat' using 1:8 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_H.dat' using 1:8 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

#set terminal postscript eps color enhanced "Arial" 30
#set output "Fr_points.eps"
#set xrange [0:10]
#set yrange [0:10]
#set ytics 0,1,10
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Froude number,  {/Arial-Italic Fr}}'
#unset title
#set key top left
#plot \
#'point_1_L.dat' using 1:8 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:8 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:8 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:8 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

set terminal postscript eps color enhanced "Arial" 30
set output "h_H_Profile-2.eps"
set xrange [0:6]
set yrange [0:0.04]
set ytics 0,0.01,0.04
#set yrange [0.0001:0.1]
#set logscale y
#set format y "10^{%L}"
#set ytics autofreq
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Thickness,  {/Arial-Italic h}_H [m]}'
unset title
set key top right
plot \
'point_1_H.dat' using 1:2 w l lt 1 lc rgb "blue" lw 6 ti "Numerical"

set terminal postscript eps color enhanced "Arial" 30
set output "h_H_Profile-3.eps"
set xrange [0:7]
set yrange [0:0.04]
#set yrange [0.0001:0.1]
#set logscale y
#set format y "10^{%L}"
#set ytics autofreq
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Thickness,  {/Arial-Italic h}_H [m]}'
unset title
set key top right
plot \
'point_2_H.dat' using 1:2 w l lt 1 lc rgb "blue" lw 6 ti "Numerical"

#set terminal postscript eps color enhanced "Arial" 30
#set output "Ri_points.eps"
#set logscale y
#set format y "10^{%L}"
#set xrange [0:10]
#set yrange [0.01:10]
#set ytics autofreq
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic Ri}}'
#unset title
#set key
#plot \
#'point_1_L.dat' using 1:9 w l lt 1 lc rgb "red" lw 6 ti "Profile 2", \
#'point_2_L.dat' using 1:9 w l lt 1 lc rgb "green" lw 6 ti "Profile 3", \
#'point_3_L.dat' using 1:9 w l lt 1 lc rgb "blue" lw 6 ti "Profile 4", \
#'point_4_L.dat' using 1:9 w l lt 1 lc rgb "black" lw 6 ti "Profile 5"

EOF

