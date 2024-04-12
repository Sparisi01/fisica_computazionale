set terminal png size 2048, 1536 font ', 36'
set output "../output/distribuzioneradiale.png"
set xrange [0:1.8]
set grid
set ylabel 'g(r)'
set xlabel 'raggio in unita di L/2'
set title "Funzione densita radiale"

plot "../data/distribuzione_radiale_gas.dat" using 1:2 w l t "Gas" lc 1, \
"../data/distribuzione_radiale_liquido.dat" using 1:2 w l t "Liquido" lc 2, \
"../data/distribuzione_radiale_solido.dat" using 1:2 w l t "Solido" lc 7
