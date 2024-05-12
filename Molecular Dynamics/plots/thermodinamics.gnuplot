set terminal png size 1600,1200 enhanced font ',28' lw 2
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set output './png/termodinamica.png'
set arrow 1 from 90, graph 0 to 90, graph 1 nohead lt 8 lw 1 dt 2

set key horizontal

set grid
set multiplot layout 2,1
set lmargin at screen 0.15
set rmargin at screen 0.95
set title "Temperatura e pressione in funzione del tempo (FCC 256)" font ",16"

#set xrange [0.4:6]
#set yrange [1.0:1.3]

set format x ""
set ylabel 'Temperatura'
plot "../data/FCC_256/thermo_gas.dat" using 1:2 w l t "ρ = 0.1", "../data/FCC_256/thermo_liquido.dat" using 1:2 w l lc 3 t "ρ = 0.8", "../data/FCC_256/thermo_solido.dat" using 1:2 w l lc 4 t "ρ = 1.2"
set ylabel 'Compressibilità'
set xlabel 'Tempo'

unset format x

unset tmargin
unset title
#set yrange [-1:6]
unset key
plot "../data/FCC_256/thermo_gas.dat" using 1:3 w l notitle, "../data/FCC_256/thermo_liquido.dat" using 1:3 w l lc 3 notitle, "../data/FCC_256/thermo_solido.dat" using 1:3 w l lc 4 notitle
