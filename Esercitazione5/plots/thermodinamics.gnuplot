set terminal pngcairo size 800,600 enhanced font 'Helvetica,14'
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set output 'termodinamica.png'
set arrow 1 from 3, graph 0 to 3, graph 1 nohead lt 8 lw 1 dt 2
set xrange [0:6]

set yrange [0.6:2]


set grid
set multiplot layout 2,1
set lmargin at screen 0.15
set rmargin at screen 0.95
set title "Temperatura e pressione in funzione del tempo" font ",16"
set tmargin at screen 0.9
#set xrange [0.4:6]
#set yrange [1.0:1.3]
set ylabel 'Temperatura'
set format x ""
plot "../data/thermo_gas.dat" using 1:2 w l t "Gas", "../data/thermo_liquido.dat" using 1:2 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:2 w l lc 4 t "Solido"
set ylabel 'Compressibilità'
set xlabel 'Tempo'
set yrange [-5:30]
unset format x
unset tmargin
unset title
#set yrange [-1:6]
unset key
plot "../data/thermo_gas.dat" using 1:3 w l t "Gas", "../data/thermo_liquido.dat" using 1:3 w l lc 3 t "Liquido", "../data/thermo_solido.dat" using 1:3 w l lc 4 t "Solido"
