set terminal png size 1600,1200 enhanced font ',28' lw 2
set output '../plots/energia.png'

set title "Evoluzione energia del sistema (FCC-500)"
set ylabel 'Energia per particella'
set xlabel '# propagazioni'

set grid
set yrange [2:-8]

set key horizontal center bottom offset 0,2

plot "../data/FCC-256/thermo_gas.dat" using ($1):($5/500) with lines t "Gas" lc 1,\
"../data/FCC-256/thermo_liquido.dat" using ($1):($5/500) with lines t "Liquido" lc 4, \
"../data/FCC-256/thermo_solido.dat" using ($1):($5/500) with lines t "Solido" lc 2

# plot "../data/FCC-256/thermo_gas.dat" using ($1):($3) with lines t "Gas" lc 1,\
# plot "../data/FCC-256/thermo_liquido.dat" using ($1):($3) with lines t "Liquido" lc 4, \
#"../data/FCC-256/thermo_solido.dat" using ($1):($3) with lines t "Solido" lc 2, \
#"../data/FCC-256-verlet/thermo_gas.dat" using ($1):($3) with lines t "Gas" lc 1, \
#"../data/FCC-256-verlet/thermo_liquido.dat" using ($1):($3) with lines t "Liquido" lc 4, \
#"../data/FCC-256-verlet/thermo_solido.dat" using ($1):($3) with lines t "Solido" lc 2
