set terminal pngcairo size 800,600 enhanced font 'Helvetica,14'
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set output 'termodinamica.png'
set arrow 1 from 3, graph 0 to 3, graph 1 nohead lt 8 lw 1 dt 2
set xrange [0:10]


set grid
set multiplot layout 3,1
set lmargin at screen 0.15
set rmargin at screen 0.95
set title "Temperatura e pressione in funzione del tempo" font ",16"

#set xrange [0.4:6]
#set yrange [1.0:1.3]
set ylabel 'Energia'
set format x ""

plot "../data/thermo_gas.dat" using 1:4 w l notitle, "../data/thermo_liquido.dat" using 1:4 w l lc 3 notitle, "../data/thermo_solido.dat" using 1:4 w l lc 4 notitle
set ylabel 'Temperatura'
unset title
plot "../data/thermo_gas.dat" using 1:2 w l t "ρ = 0.1", "../data/thermo_liquido.dat" using 1:2 w l lc 3 t "ρ = 0.8", "../data/thermo_solido.dat" using 1:2 w l lc 4 t "ρ = 1.2"
set ylabel 'Compressibilità'
set xlabel 'Tempo'

unset format x
set xtics 1

unset tmargin
unset title
#set yrange [-1:6]
unset key
plot "../data/thermo_gas.dat" using 1:3 w l notitle, "../data/thermo_liquido.dat" using 1:3 w l lc 3 notitle, "../data/thermo_solido.dat" using 1:3 w l lc 4 notitle
