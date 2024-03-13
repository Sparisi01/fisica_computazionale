reset
set datafile commentschars '#@&'
set terminal png size 2048, 1536 font ', 36'
set output 'plot_errore_derivata.png'
set key top right
set grid
set title 'Errore in funzione del passo h' font ', 50'
set xrange [1e-17:10]
set yrange [1e-8:80]
set logscale x
set logscale y
set format y '1x10^{%T}'
set xlabel 'h'
set ylabel 'Errore'
e = exp(1)
error_double = 1.11e-16
error_float = 5.96e-08
f(x) = e*(error_double/x + x/2)
g(x) = e*(error_float/x + x/2)
plot 'dati_double.dat' using 1 : 2 with linespoints pointtype 6 pointsize 4 linewidth 3.5 linecolor 'red' title 'Doppia precisione', f(x) notitle with lines linewidth 3.5 linecolor 'red', 'dati_float.dat' using 1 : 2 with linespoints pointtype 6 pointsize 4 linewidth 3.5 linecolor 'blue' title 'Singola precisione', g(x) notitle with lines linewidth 3.5 linecolor 'blue'
