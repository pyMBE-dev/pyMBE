set terminal png size 1280,820
set term png font "helvetica,26"

set output 'average_charge.png'

set border 3
unset tics
set xlabel 'pH'
set ylabel 'Average charge'
set xtics border
set xrange [0:14]
set ytics axis
set tics out

a(x) = 0

plot "prot_charge.dat" u 1:2 w lp pt 5 ps 2 lw 2 lc rgb "blue" t "Histadin-5",\
"data.dat" u 1:2 w lp pt 5 ps 2 lw 2 lc rgb "green" t "Original data",\
"prot_charge.dat" u 1:2:3 w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" notitle,\
[0:14] a(x) w l lw 2 lc rgb "black" lt '-' notitle

