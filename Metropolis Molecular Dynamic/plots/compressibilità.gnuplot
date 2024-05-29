set terminal png size 1600,1200 enhanced font ',28' lw 2
set output '../plots/compressibilita.png'

set title "Evoluzione energia del sistema (FCC-500)"
set ylabel 'Energia per particella'
set xlabel '# propagazioni'

set grid

set key horizontal center bottom offset 0,2


plot "../data/FCC-256/thermo_liquido.dat" using 1:4 with lines t "Liquido" lc 4, \
"../data/FCC-256-verlet/thermo_liquido.dat" using ($1/1e-3):4 with lines t "Liquido" lc 5, \
