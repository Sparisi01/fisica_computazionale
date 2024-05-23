#set terminal png size 2048, 1536 font ', 36'
#set output "../output/distribuzioneradiale.png"

# Set terminal to PNG or PDF for high-quality output
set terminal png size 1600,1200 enhanced font ',28' lw 2
# For PDF output, use:
# set terminal pdfcairo size 8cm,6cm font 'Arial,10'

# Set output file name
set output './png/peso_boltzmann.png'
# For PDF output, change the file extension to '.pdf'
set ylabel "Densità radiale g(r)"

set lmargin at screen 0.15
set rmargin at screen 0.95


set tmargin at screen 0.9
#set xrange [0:1]
set grid

set ytics 0.5
set format x ""
set title "Confronto funzione densità radiale e peso di Boltzmann" font ",32"

set style line 1 lc rgb '#0060ad' lt 1 lw 4 pt 7 ps 0.5
set style line 2 lc rgb '#dd181f' lt 1 lw 4 pt 5 ps 0.5

set xrange [0:1]

set xlabel 'r'
set key box
set key spacing 2 offset -1,-1

f(x) = exp(-((4*((1/(x*29.47/2))**12 - (1/(x*29.47/2))**6)))/1.16)

plot "../data/FCC_256/distribuzione_radiale_gas.dat" using 2:3 with l t "FCC 256" lc 1, \
f(x) t "Peso di Boltzmann" lw 2

