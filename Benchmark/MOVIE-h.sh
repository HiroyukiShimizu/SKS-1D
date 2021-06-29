#! /bin/bash

horizontal_axis="x"
#horizontal_axis="r"

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

#set size 1.823, 1.5
#set origin 0.177,1.55
set size 1.75, 1.5
set origin 0.25,1.55
set xrange [0:30]
set xtics 0,5,30
set format x ""
set yrange [0:5]
set ytics 0,1,5
unset xlabel #'{/=45 Distance,  {/Arial-Italic ${horizontal_axis}} [m]}'
#set ylabel '{/=45 {/Arial-Italic h} [m]}'
set ylabel '{/=45 Thickness [m]}'
set title '{/=45 Time,  {/Arial-Italic t} = ${t} [s]}'
set arrow 1 from 9.68,0 to 9.68,5 nohead lt 0 lc rgb "gray" lw 10
set arrow 2 from 2.65,0 to 2.65,5 nohead lt 0 lc rgb "gray" lw 10
set arrow 3 from 7.78,0 to 7.78,5 nohead lt 0 lc rgb "gray" lw 10
#set arrow 4 from 13.68,0 to 13.68,5 nohead lt 0 lc rgb "gray" lw 10
set key
plot \
'prim.L${number}' using 20:2 w l lt 1 lc rgb "red" lw 6 ti "Dilute current"

set size 1.862, 0.7
set origin 0.138,0.9
set xrange [0:30]
set xtics 0,5,30
set format x ""
set yrange [0:0.04]
set ytics 0,0.02, 0.04
unset xlabel #'{/=45 Distance,  {/Arial-Italic ${horizontal_axis}} [m]}'
#set ylabel '{/=45 {/Arial-Italic h}_H [m]}'
set ylabel '{/=45 Thcikness [m]}'
unset title
set arrow 1 from 9.68,0 to 9.68,0.04 nohead lt 0 lc rgb "gray" lw 10
set arrow 2 from 2.65,0 to 2.65,0.04 nohead lt 0 lc rgb "gray" lw 10
set arrow 3 from 7.78,0 to 7.78,0.04 nohead lt 0 lc rgb "gray" lw 10
#set arrow 4 from 13.68,0 to 13.68,0.04 nohead lt 0 lc rgb "gray" lw 10
set key
plot \
'prim.H${number}' using 8:2 w l lt 1 lc rgb "blue" lw 6 ti "Basal current (Bedload)"

set size 1.9, 0.8
set origin 0.1,0.15
set xrange [0:30]
set xtics 0,5,30
set format x "%g"
set yrange [0:0.004]
set ytics 0,0.002, 0.004
set xlabel '{/=45 Distance,  {/Arial-Italic ${horizontal_axis}} [m]}'
#set ylabel '{/=45 {/Arial-Italic z}_b [m]}'
set ylabel '{/=45 Thickness [m]}'
unset title
set arrow 1 from 9.68,0 to 9.68,0.004 nohead lt 0 lc rgb "gray" lw 10
set arrow 2 from 2.65,0 to 2.65,0.004 nohead lt 0 lc rgb "gray" lw 10
set arrow 3 from 7.78,0 to 7.78,0.004 nohead lt 0 lc rgb "gray" lw 10
#set arrow 4 from 13.68,0 to 13.68,0.004 nohead lt 0 lc rgb "gray" lw 10
#set label 1 at graph 0.02,5.88 "{/=30 Profile 2}"
#set label 2 at graph 0.18,5.88 "{/=30 Profile 3}"
#set label 3 at graph 0.31,5.88 "{/=30 Break of slope}"
#set label 4 at graph 0.45,5.88 "{/=30 Confine/Unconfined boundary}"
set key
plot \
'prim.D${number}' using 4:2 w l lt 1 lc rgb "black" lw 6 ti "Deposit"



unset multiplot


set out
EOF

done

convert -delay 5 -background white -alpha off h_????.eps h.gif

