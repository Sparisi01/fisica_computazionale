#set terminal png size 2048, 1536 font ', 36'
#set output "../output/distribuzioneradiale.png"

set lmargin at screen 0.15
set rmargin at screen 0.95
set multiplot layout 3,1
set xrange [0:1.8]
set grid
set xtics 0.1
set ytics 0.5
set format x ""
set title "Funzione densita radiale"


plot "../data/distribuzione_radiale_gas.dat" using 1:2 with boxes t "Gas" lc 1
unset title
plot "../data/distribuzione_radiale_liquido.dat" using 1:2 with boxes t "Liquido" lc 2
set xlabel 'raggio in unita di L/2'
unset format x
set ytics 1
plot "../data/distribuzione_radiale_solido.dat" using 1:2 with boxes t "Solido" lc 7
