#set terminal png size 2048, 1536 font ', 36'
#set output "../output/distribuzioneradiale.png"

# Set terminal to PNG or PDF for high-quality output
set terminal pngcairo size 800,600 enhanced font 'Helvetica,14'
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set output './png/distribuzione_radiale.png'
# For PDF output, change the file extension to '.pdf'
set ylabel "Densità radiale g(r)" offset -2,-9.6

set lmargin at screen 0.15
set rmargin at screen 0.95

set multiplot layout 3,1
set tmargin at screen 0.9
#set xrange [0:1]
set grid

set ytics 0.5
set format x ""
set title "Funzione densità radiale" font ",16"

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 0.5

set xrange [0:1]


f(x) = exp(-((4*((1/(x*29.47/2))**12 - (1/(x*29.47/2))**6))-0.01)/1.16)

set y2label "GAS"
plot "../data/FCC_256/distribuzione_radiale_gas.dat" using 2:3 with l t "FCC 256" lc 1, \
"../data/BCC_128/distribuzione_radiale_gas.dat" using 2:3 with l t "BCC 128" lc 1, \
f(x) t "Peso di Boltzman"
unset ylabel
unset tmargin
set y2label "LIQUIDO"
unset title
set yrange [0:3]


plot "../data/FCC_256/distribuzione_radiale_liquido.dat" using 2:3 with l t "FCC 256" lc 2, \
"../data/BCC_128/distribuzione_radiale_liquido.dat" using 2:3 with l t "BCC 128" lc 2
set xlabel 'Raggio in unità di L/2'
unset format x
set ytics 1
set y2label "SOLIDO"
set yrange [0:6]

plot "../data/FCC_256/distribuzione_radiale_solido.dat" using 2:3 with l t "FCC 256" lc 7, \
"../data/BCC_128/distribuzione_radiale_solido.dat" using 2:3 with l t "BCC 128" lc 7

