# Set terminal to PNG or PDF for high-quality output
set terminal png size 1600,1200 enhanced font ',28' lw 4
set output './compressibilita.png'

set yrange [-1:18]

set tmargin at screen 0.9
# Set title and labels
set grid
set lmargin at screen 0.15
set rmargin at screen 0.95
set xrange [-0.1:1.3]

set xlabel 'Densità'
set ylabel 'Compressibilità'
set title "Compressibilità in funzione della densità FCC-256)" font ",32"

# Set plot style and colors
# Add more styles if needed
# Set grid
set grid
set key left spacing 2 offset 1,-1
set key box
# Plot data
plot "../data/FCC-256/pressione_temperatura_rho.dat" using 1:3 w p lc 1 ps 4 pt 6 t "metropolis", \
"../../Molecular Dynamics/data/FCC_256/pressione_temperatura_rho.dat" using 1:3 w p lc 4 ps 4 pt 8 t "verlet"

# Add legend
set key top left

# Add any other customizations
