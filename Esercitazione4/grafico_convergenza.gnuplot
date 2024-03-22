reset
set monochrome
set ytics nomirror
set xtics nomirror
set key right
set key box lt -1 lw 1
set key spacing 2 font "Helvetica, 12"

set grid

set xlabel 'passo h' 
set ylabel 'Differenza passo precedente' font "Helvetica, 12"
set title "Convergenza al variare del passo h" font "Helvetica, 14"

set logscale y
set logscale x
set format y '10^{%T}'
set format x '10^{%T}'

set xrange [1E-1:1E-8];
set yrange [1E-7:1];

plot "convergenza_eulero.dat" u 1:2 pt 6 t "Metodo Eulero", \
"convergenza_rk4.dat" u 1:2 pt 4 t "Metodo ek4" 