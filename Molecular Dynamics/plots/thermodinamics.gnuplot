
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set terminal png size 1600,1200 enhanced font ',28' lw 2
set output './png/termostato.png'

set arrow 1 from 90, graph 0 to 90, graph 1 nohead lt 8 lw 1 dt 2

set yrange [0.9:1.4]

set key horizontal

set grid
set multiplot layout 2,1
set lmargin at screen 0.15
set rmargin at screen 0.95
set title "Fluttuazioni in energia con termostato (FCC 256)" font ",32"

#set xrange [0.4:6]
#set yrange [1.0:1.3]

set format x ""
set ylabel 'Temperatura'
plot "../data/thermo_gas.dat" using 1:2 w l t "ρ = 0.1", "../data/thermo_liquido.dat" using 1:2 w l lc 3 t "ρ = 0.8", "../data/thermo_solido.dat" using 1:2 w l lc 4 t "ρ = 1.2"
set ylabel 'Energia per particella'
set xlabel 'Tempo'
unset yrange
unset format x
unset tmargin
unset title
#set yrange [-1:6]
unset key

plot "../data/thermo_gas.dat" using 1:($5/256) w l notitle, "../data/thermo_liquido.dat" using 1:($5/256) w l lc 3 notitle, "../data/thermo_solido.dat" using 1:($5/256) w l lc 4 notitle
