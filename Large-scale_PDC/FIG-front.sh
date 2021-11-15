#! /bin/bash

#horizontal_axis="x"
horizontal_axis="r"

gnuplot << EOF
set terminal postscript eps color enhanced "Arial" 30
set output "Front-position.eps"
set xrange [0:200]
set yrange [0:10000]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 Front position,  {/Arial-Italic ${horizontal_axis}}  [m]}'
unset title
set key bottom right
plot \
'front_H.dat' using 1:2 w l lt 1 lc rgb "blue" lw 6 ti "Dense current", \
'front_L.dat' using 1:2 w l lt 1 lc rgb "red" lw 6 ti "Dilute current"

set terminal postscript eps color enhanced "Arial" 30
set output "h_FC.eps"
set xrange [0:200]
set yrange [0:800]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic h}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:3 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic h}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "rho_FC.eps"
set xrange [0:200]
set yrange [0:5]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Symbol-Oblique r}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:4 w l lt 1 lc rgb "red" lw 6 ti "{/Symbol-Oblique r}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "u_FC.eps"
set xrange [0:200]
set yrange [0:100]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic u}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:5 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic u}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "n_aFC.eps"
set xrange [0:200]
set yrange [0:1]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic n}_{a,FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:6 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic n}_{a,FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "n_sFC.eps"
set xrange [0:200]
set yrange [0:1]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic n}_{s,FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:7 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic n}_{s,FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "T_FC.eps"
set xrange [0:200]
set yrange [0:1000]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic T}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:8 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic T}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "rhoFC_rhoa.eps"
set xrange [0:200]
set yrange [0:4]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Symbol-Oblique r}_{FC}/{/Symbol-Oblique r}_a({/Arial-Italic t})}'
unset title
unset key
plot \
'front_L.dat' using 1:11 w l lt 1 lc rgb "red" lw 6 ti "{/Symbol-Oblique r}_{FC}/{/Symbol-Oblique r}_a({/Arial-Italic t})", \
1 w l lt 0 lc rgb "black" lw 6 ti "1"

set terminal postscript eps color enhanced "Arial" 30
set output "Fr_FC.eps"
set xrange [0:200]
set yrange [0:5]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic Fr}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:9 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic Fr}_{FC}({/Arial-Italic t})"

set terminal postscript eps color enhanced "Arial" 30
set output "Ri_FC.eps"
set logscale y
set format y "10^{%L}"
set xrange [0:200]
set yrange [0.001:1000000]
set xlabel '{/=35 Time,  {/Arial-Italic t}  [s]}'
set ylabel '{/=35 {/Arial-Italic Ri}_{FC}({/Arial-Italic t})}'
unset title
unset key
plot 'front_L.dat' using 1:10 w l lt 1 lc rgb "red" lw 6 ti "{/Arial-Italic Ri}_{FC}({/Arial-Italic t})"

EOF

