reset
set grid
set xlabel "Raggio [R0]"
set ylabel "Massa [M0]"
set logscale x
set logscale y
plot "dati_stelle1.dat" using 1:2 title "Stelle tipo 1", \
"dati_stelle2.dat" using 1:2 title "Stelle tipo 2", \
"dati_stelle3.dat" using 1:2 title "Stelle tipo 3"