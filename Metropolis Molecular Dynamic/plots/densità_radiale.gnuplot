set multiplot layout 3,1

set xrange [0:1]

f(x) = exp(-((4*((1/(x*29.47/2))**12 - (1/(x*29.47/2))**6)))/1.1)
plot "../data/distribuzione_radiale_gas.dat" using 2:3 w lines t "Gas", f(x)


plot "../data/distribuzione_radiale_liquido.dat" using 2:3 w lines t "Liquido"

plot "../data/distribuzione_radiale_solido.dat" using 2:3 w lines t "Solido"