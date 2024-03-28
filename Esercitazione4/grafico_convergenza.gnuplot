reset

#set terminal tikz size 3.5in,2.4in monochrome standalone font "Helvetica,10"
set monochrome
set datafile commentschars '#@&'
#set output 'convergenza_algoritmi.tex'


set ytics nomirror
set xtics nomirror
set key right
set key box lt -1 lw 1
set key spacing 2 font ", 12"

set grid

set xlabel 'passo h' font ", 12"
set ylabel 'Differenza passo precedente' font ", 12"
set title "Convergenza al variare del passo h" font ", 14"

set logscale y
set logscale x
set format y '10^{%T}'
set format x '10^{%T}'

set xrange [1E-1:1E-8];
set yrange [1E-7:1];


plot "./data/convergenza_eulero.dat" u 1:2 pt 7 t "Metodo Eulero", \
"./data/convergenza_rk4.dat" u 1:2 pt 4 t "Metodo ek4" 
