reset
set datafile commentschars '#@&'
set terminal png size 2048, 1536 font ', 36'
set output 'plot_psi.png'
set grid
set title "Modulo quadro funzione d'onda in funzione di x" font ', 50'
set xrange [0:3]
set yrange [0:1.5]
set xlabel 'x'
set ylabel 'Densità probabilità'
set arrow from 0.573962, graph 0 to 0.573962, graph 1 nohead lt 1 lc rgb "red" lw 3
set arrow from (0.573962-0.518406), graph 0 to (0.573962-0.518406), graph 1 nohead lw 2 lt 0 lc rgb "red"
set arrow from (0.573962+0.518406), graph 0 to (0.573962+0.518406), graph 1 nohead lw 2 lt 0 lc rgb "red"
plot "valoriPsi.dat" with lines lw 3 notitle