



# Set terminal to PNG or PDF for high-quality output
set terminal pngcairo size 800,600 enhanced font 'Helvetica,12'
set output 'compressibilita.png'

set tmargin at screen 0.9
# Set title and labels
set grid
set lmargin at screen 0.15
set rmargin at screen 0.95
set xrange [-0.1:1.2]
set yrange [0:3]
set xlabel 'Densità' font ",14"
set ylabel 'Compressibilità' font ",14"
set title "Compressibilità in funzione della densità" font ",16"

# Set plot style and colors
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb '#dd181f' lt 1 lw 2 pt 5 ps 0.5
# Add more styles if needed
# Set grid
set grid



# Plot data
plot "../data/pressione_temperatura_rho.dat" using 1:3 w lp ps 1 pt 5 notitle
# Add more data sets as needed

# Add legend
set key top left

# Add any other customizations
