set grid
set multiplot layout 2,1
set lmargin at screen 0.15
set rmargin at screen 0.95
set ylabel 'Temperatura'
set format x ""
plot "../data/thermo_gas.dat" using 1:2 w l t "Gas", "../data/thermo_liquido.dat" using 1:2 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:2 w l lc 4 t "Solido"
set ylabel 'Pressione'
set xlabel 'Tempo'

unset format x
plot "../data/thermo_gas.dat" using 1:3 w l t "Gas", "../data/thermo_liquido.dat" using 1:3 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:3 w l lc 4 t "Solido"
