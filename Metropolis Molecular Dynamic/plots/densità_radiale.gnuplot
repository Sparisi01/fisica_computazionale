set terminal png size 1600,1200 enhanced font ',28' lw 4
set output '../plots/densita_radiale.png'

set multiplot layout 3,1

set xrange [0:1]
set xtics 0.1
set mxtics 2

set lmargin at screen 0.15
set rmargin at screen 0.95

set title "Distribuzione radiale (FCC-256)"
set format x ""
L = 36.840315
#f(x) = exp(-((4*((1/(x*L/2))**12 - (1/(x*L/2))**6)))/1.1) 
plot "../data/FCC-256/distribuzione_radiale_gas.dat" using 2:3 w lines t "Metropolis" lc 1, \
"../data/FCC-256-verlet/distribuzione_radiale_gas.dat" using 2:3 w p t "Verlet" lc 4
unset title

plot "../data/FCC-256/distribuzione_radiale_liquido.dat" using 2:3 w lines notitle lc 1, \
"../data/FCC-256-verlet/distribuzione_radiale_liquido.dat" using 2:3 notitle lc 4

unset format x
set ylabel "Densità radiale g(r)" offset -2,9.6
set xlabel 'Raggio in unità di L/2'
plot "../data/FCC-256/distribuzione_radiale_solido.dat" using 2:3 w lines notitle lc 1, \
"../data/FCC-256-verlet/distribuzione_radiale_solido.dat" using 2:3 notitle lc 4