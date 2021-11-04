set terminal png size 1280,820
set term png font "helvetica,26"

set output 'average_rg.png'

set border 3
unset tics
set xlabel 'pH'
set ylabel 'Average radius of gyration / nm'
set xtics border
set xrange [0:14]
set ytics axis
set tics out

a(x) = 0

plot "prot_rg.dat" u 1:2 w lp pt 5 ps 2 lw 2 lc rgb "blue" t "Histadin-5",\
"prot_rg.dat" u 1:2:3 w yerrorbars pt 5 ps 2 lw 2 lc rgb "blue" notitle,\
"data.dat" u 1:($5/10):6 w yerrorbars pt 5 ps 2 lw 2 lc rgb "green" t "Original data",\
"data.dat" u 1:($5/10) w l lw 2 lc rgb "green" notitle

