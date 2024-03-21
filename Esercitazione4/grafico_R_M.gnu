reset
set multiplot layout 1,2
set grid
set xlabel "Raggio [R0]"
set ylabel "Massa [M0]"
plot "dati_stelle1.dat" u 1:2 w lp t "Stelle tipo 1 mr4", \
"dati_stelle1.dat" u 3:4 w lp t "Stelle tipo 1 eulero"

set xlabel "Delta Raggio [R0]"
set ylabel "Delta Massa [M0]"

plot "dati_stelle1.dat" u 5:6

