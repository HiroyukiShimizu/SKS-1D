#! /bin/bash

horizontal_axis="x"
#horizontal_axis="r"

files=`ls prim.L????`

for file in ${files};do 
t=`grep   "###" ${file} |awk '{printf $3}'`
number=${file##*L}
gnuplot << EOF

set terminal postscript eps color enhanced "Arial" 30
set output "deposit_mass_${number}.eps"
#set size 2, 1
set size 1.5, 1
set xrange [0:13.68]
set xtics 0,5,15
set yrange [0:6]
#set yrange [0.1:100]
#set logscale y
#set format y "10^{%L}"
#set ytics autofreq
set xlabel '{/=35 Distance,  {/Arial-Italic ${horizontal_axis}} [m]}'
set ylabel '{/=35 Deposit mass,  {/Symbol-Oblique f}_{sD}{/Symbol-Oblique r}_s{/Arial-Italic z}_b [kg/m^2]}'
set title '{/=35 Time,  {/Arial-Italic t} = ${t} [s]}'
set arrow 1 from 9.68,0 to 9.68,6 nohead lt 0 lc rgb "gray" lw 10
set arrow 2 from 2.65,0 to 2.65,6 nohead lt 0 lc rgb "gray" lw 10
set arrow 3 from 7.78,0 to 7.78,6 nohead lt 0 lc rgb "gray" lw 10
#set arrow 4 from 13.68,0 to 13.68,6 nohead lt 0 lc rgb "gray" lw 10
set key
plot \
'prim.D${number}' using 4:5 w l lt 1 lc rgb "black" lw 6 ti "Numerical"

set out
EOF

done

convert -delay 5 -background white -alpha off deposit_mass_????.eps deposit_mass.gif

