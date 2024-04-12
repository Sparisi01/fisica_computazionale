reset

#set terminal tikz size 3.5in,2.4in monochrome standalone font "Helvetica,10"
set terminal png size 2048, 1536 font ', 36'
set monochrome
set datafile commentschars '#@&'
#set output 'convergenza_algoritmi.tex'
set output "./out/grafico_convergenza.png"


set ytics nomirror
set xtics nomirror
set key right
set key box lt -1 lw 1
set key spacing 2 font ", 36"

set grid

set xlabel 'passo h' font ", 36"
set ylabel 'Differenza passo h/2' font ", 36"
set title "Convergenza al variare del passo h" font ", 36"

set logscale y
set logscale x
set format y '10^{%T}'
set format x '10^{%T}'

set xrange [1E-1:1E-8];
set yrange [1E-7:1];


plot "./data/convergenza_eulero.dat" u 1:2 pt 7 pointsize 4 linewidth 3.5 t "Metodo Eulero", \
"./data/convergenza_rk4.dat" u 1:2 pt 4 pointsize 4 linewidth 3.5 t "Metodo rk4" 


