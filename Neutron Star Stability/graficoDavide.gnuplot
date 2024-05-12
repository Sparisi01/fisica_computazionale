



# Set terminal to PNG or PDF for high-quality output
set terminal pngcairo size 800,600 enhanced font 'Helvetica,14'
set output 'graficoerrorerkeulero.png'

set tmargin at screen 0.9
# Set title and labels
set grid
set lmargin at screen 0.15
set rmargin at screen 0.95


set xlabel 'h/h₀'
set ylabel 'Errore rispetto a raggio di riferimento'
set title "Errore metodo di Eulero e Runge kutta al variare del passo h" 


set grid

set logscale y
set logscale x



# Plot data
plot "./data/errore_rk.dat" using 1:3 w lp ps 1 pt 5 t "Runge Kutta", \
"./data/errore_eu.dat" using 1:3 w lp ps 1 pt 5 t "Eulero",
# Add more data sets as needed
# Add legend
set key top right spacing 2

# Add any other customizations
