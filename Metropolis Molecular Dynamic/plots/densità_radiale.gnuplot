set terminal png size 1600,1200 enhanced font ',28' lw 2
set output '../plots/densita_radiale.png'

set multiplot layout 3,1

set xrange [0:1]

set lmargin at screen 0.15
set rmargin at screen 0.95

set title "Distribuzione radiale"
set format x ""
L = 36.840315
f(x) = exp(-((4*((1/(x*L/2))**12 - (1/(x*L/2))**6)))/1.1) 
plot "../data/distribuzione_radiale_gas.dat" using 2:3 w lines t "Gas" lc 1 lw 2, f(x) w l 

unset title

plot "../data/distribuzione_radiale_liquido.dat" using 2:3 w lines t "Liquido" lc 2 lw 2

unset format x
set ylabel "Densità radiale g(r)" offset -2,9.6
set xlabel 'Raggio in unità di L/2'
plot "../data/distribuzione_radiale_solido.dat" using 2:3 w lines t "Solido" lc 4 lw 2