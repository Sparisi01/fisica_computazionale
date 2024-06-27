set terminal png size 1600,1200 enhanced font ',28' lw 4
set output '../plots/compressibilita.png'
set datafile separator ";"
set title "Evoluzione comprssibilità del sistema (FCC-256)"
set ylabel 'Compressibilità'
set xlabel '# propagazioni'

set grid

set key horizontal center bottom offset 0,2

plot "../data/FCC-256/thermo_01.csv" using ($1):($4) with lines t "Gas" lc 1 lw 0.5,\
"../data/FCC-256/thermo_02.csv" using ($1):($4) with lines t "Liquido" lc 4 lw 0.5, \
"../data/FCC-256/thermo_03.csv" using ($1):($4) with lines t "Solido" lc 2 lw 0.5


