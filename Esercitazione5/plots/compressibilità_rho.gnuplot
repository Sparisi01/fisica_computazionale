set grid
set xrange [-0.1:1.5]
set yrange [0:2]
set xlabel 'ρ'
set ylabel 'P'

plot "../data/pressione_temperatura_rho" using 1:3 ps 1 pt 5 notitle
