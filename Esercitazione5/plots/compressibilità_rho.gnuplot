set grid
set xrange [-0.1:1.5]
set yrange [0:3]
set xlabel 'Densita'
set ylabel 'Compressibilita'

plot "../data/pressione_temperatura_rho.dat" using 1:3 ps 1 pt 5 notitle
