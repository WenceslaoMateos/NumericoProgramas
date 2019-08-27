set autoscale
unset log
unset label
set xtic auto
set ytic auto

set title "Euler"
set xlabel "x"
set ylabel "f(x)"

plot "datos.dat" using 1:2 title 'euler' with linespoints



