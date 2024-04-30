set grid
set multiplot layout 2,1
set lmargin at screen 0.15
set rmargin at screen 0.95
#set xrange [0.4:6]
#set yrange [1.0:1.3]
set ylabel 'Temperatura'
set format x ""
plot "../data/thermo_gas.dat" using 1:2 w l t "Gas", "../data/thermo_liquido.dat" using 1:2 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:2 w l lc 4 t "Solido"
set ylabel 'Pressione'
set xlabel 'Tempo'

unset format x
#set yrange [-1:6]

plot "../data/thermo_gas.dat" using 1:3 w l t "Gas", "../data/thermo_liquido.dat" using 1:3 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:3 w l lc 4 t "Solido"
