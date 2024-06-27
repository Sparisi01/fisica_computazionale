set terminal png size 1600,1200 enhanced font ',28' lw 4
set output '../plots/energia.png'
set datafile separator ";"
set title "Evoluzione energia del sistema (FCC-256)"
set ylabel 'Energia per particella'
set xlabel '# propagazioni'

set grid
set yrange [2:-8]

set key horizontal center bottom offset 0,2

plot "../data/FCC-256/thermo_01.csv" using ($1):($5/256) with lines t "Gas" lc 1 lw 0.5,\
"../data/FCC-256/thermo_02.csv" using ($1):($5/256) with lines t "Liquido" lc 4 lw 0.5, \
"../data/FCC-256/thermo_03.csv" using ($1):($5/256) with lines t "Solido" lc 2 lw 0.5

