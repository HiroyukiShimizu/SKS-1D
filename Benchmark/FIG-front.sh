#! /bin/bash

horizontal_axis="x"
#horizontal_axis="r"

gnuplot << EOF
set terminal postscript eps color enhanced "Arial" 30
set output "Front-position.eps"
set xrange [0:7]
#set xrange [0:4]
#set xtics 0,2,10
#set xtics 0,1,5
set yrange [0:20]
#set yrange [0:13.68]
#set ytics 0,5,30
#set ytics 0,5,15
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Front position,  {/Arial-Italic ${horizontal_axis}}_N [m]}'
unset title
set key bottom right
plot \
9.68 w l lt 0 lc rgb "gray" lw 6 notitle, \
2.65 w l lt 0 lc rgb "gray" lw 6 notitle, \
7.78 w l lt 0 lc rgb "gray" lw 6 notitle, \
'front_L.dat' using 1:14 w l lt 1 lc rgb "red" lw 6 ti "Dilute current", \
'front_H.dat' using 1:4 w l lt 1 lc rgb "blue" lw 6 ti "Basal current (Bedload)", \
'front_H.dat' using 1:5 w l lt 0 lc rgb "black" lw 6 ti "Deposit"


#set terminal postscript eps color enhanced "Arial" 30
#set output "h_FC.eps"
#set xrange [0:10]
#set yrange [0:3]
#set ytics 0,1,3
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic h}_{FC} [m]}'
#unset title
#unset key
#plot 'front_L.dat' using 1:3 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic h}_{FC}({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "rho_FC.eps"
#set xrange [0:10]
#set yrange [0:3]
#set ytics 0,0.5,3
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Symbol-Oblique r}_{FC} [kg/m^3]}'
#unset title
#unset key
#plot 'front_L.dat' using 1:4 w l lt 1 lc rgb "red" lw 6 ti "{/Symbol-Oblique r}_{FC}({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "u_FC.eps"
#set xrange [0:5]
#set xtics 0,1,5
#set yrange [0:8]
#set ytics 0,2,8
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 Frontal velocity,  {/Arial-Italic u}_{N} [m/s]}'
#unset title
#unset key
#plot \
#'front_L.dat' using 1:5 w l lt 1 lc rgb "black" lw 6 ti "{/Arial-Italic u}_N({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "n_aFC.eps"
#set xrange [0:10]
#set yrange [0:1]
#set ytics 0,0.1,1
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic n}_{a,FC}}'
#unset title
#unset key
#plot 'front_L.dat' using 1:6 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic n}_{a,FC}({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "n_sFC.eps"
#set xrange [0:10]
#set yrange [0:1]
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic n}_{s,FC}}'
#unset title
#unset key
#plot 'front_L.dat' using 1:7 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic n}_{s,FC}({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "T_FC.eps"
#set xrange [0:10]
#set yrange [0:400]
#set ytics 0,50,400
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic T}_{FC} [K]}'
#unset title
#set key
#plot \
#284 w l lt 0 lc rgb "black" lw 6 ti "{/Arial-Italic T}_{a}", \
#'front_L.dat' using 1:8 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic T}_{FC}"

#set terminal postscript eps color enhanced "Arial" 30
#set output "rhoFC_rhoa.eps"
#set xrange [0:10]
#set yrange [1:2]
#set ytics 0,0.1,2
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Symbol-Oblique r}_{FC}/{/Symbol-Oblique r}_a}'
#unset title
#unset key
#plot 'front_L.dat' using 1:11 w l lt 1 lc rgb "red" lw 6 ti "{/Symbol-Oblique r}_{FC}/{/Symbol-Oblique r}_a({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "Fr_FC.eps"
#set xrange [0:10]
#set yrange [0:4]
#set ytics 0,0.5,4
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic Fr}_{FC}}'
#unset title
#unset key
#plot 'front_L.dat' using 1:9 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic Fr}_{FC}({/Arial-Italic t})"

#set terminal postscript eps color enhanced "Arial" 30
#set output "Ri_FC.eps"
#set logscale y
#set format y "10^{%L}"
#set xrange [0:10]
#set yrange [0.01:10]
#set ytics autofreq
#set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
#set ylabel '{/=35 {/Arial-Italic Ri}_{FC}}'
#unset title
#unset key
#plot 'front_L.dat' using 1:10 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic Ri}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "Log_front-position.eps"
set logscale x
set format x "10^{%L}"
set logscale y
set format y "10^{%L}"
set xrange [0.1:10]
set xtics autofreq
set yrange [0.1:20]
set ytics autofreq
set xlabel '{/=35 Time,  {/Arial-Italic t} [s]}'
set ylabel '{/=35 Front position,  {/Arial-Italic ${horizontal_axis}} [m]}'
unset title
set key bottom right
plot \
9.68 w l lt 0 lc rgb "gray" lw 6 notitle, \
2.65 w l lt 0 lc rgb "gray" lw 6 notitle, \
7.78 w l lt 0 lc rgb "gray" lw 6 notitle, \
'front_L.dat' using 1:14 w l lt 1 lc rgb "red" lw 6 ti "Dilute current", \
'front_H.dat' using 1:4 w l lt 1 lc rgb "blue" lw 6 ti "Basal current (Bedload)", \
'front_H.dat' using 1:5 w l lt 0 lc rgb "black" lw 6 ti "Deposit"


EOF

