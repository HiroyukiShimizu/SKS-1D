#! /bin/bash

#horizontal_axis="x"
horizontal_axis="r"

files=`ls prim.L????`

for file in ${files};do 
t=`grep   "###" ${file} |awk '{printf $3}'`
number=${file##*L}
gnuplot << EOF

set terminal postscript eps color enhanced "Arial" 45
set output "h_${number}.eps"
set size 2, 3
set origin 0.1,0.1
set multiplot layout 3,1

set size 1.9, 1.5
set origin 0.1,1.55
set xrange [0:10000]
set xtics 0,2000,10000
set format x ""
set yrange [0:800]
set ytics 0,200,800
unset xlabel
set ylabel '{/=45 {/Arial-Italic h}  [m]}'
set title '{/=45 Time,  {/Arial-Italic t} = ${t}  [s]}'
set key
plot \
'prim.L${number}' using 1:2 w l lt 1 lc rgb "red" lw 6 ti "Dilute current"

set size 1.9, 0.6
set origin 0.1,1.0
set xrange [0:10000]
set xtics 0,2000,10000
set format x ""
set yrange [0:0.5]
set ytics 0,0.5,0.5
unset xlabel #'{/=45 Distance,  {/Arial-Italic ${horizontal_axis}}  [m]}'
set ylabel '{/=45 {/Arial-Italic h}_H  [m]}'
unset title
set key
plot \
'prim.H${number}' using 1:2 w l lt 1 lc rgb "blue" lw 6 ti "Dense current"

set size 1.9, 0.9
set origin 0.1,0.15
set xrange [0:10000]
set xtics 0,2000,10000
set format x "%g"
set yrange [0:1]
set ytics 0,0.5,1
set xlabel '{/=45 Distance,  {/Arial-Italic ${horizontal_axis}}  [m]}'
set ylabel '{/=45 {/Arial-Italic z}_b  [m]}'
unset title
set key
plot \
'prim.D${number}' using 1:2 w l lt 1 lc rgb "black" lw 6 ti "Deposit"

unset multiplot

set out
EOF

done

convert -delay 5 -background white -alpha off h_????.eps h.gif

