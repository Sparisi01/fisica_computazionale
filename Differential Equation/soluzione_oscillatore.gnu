reset
set multiplot layout 1,2

set grid

plot "plot_E_esplicito.dat" using 1:2 with lines title "Eulero Esplicito", \
"plot_E_implicito.dat" using 1:2 with lines title "Eulero Implicito", \
"plot_Trapezio.dat" using 1:2 with lines title "Trapezi", \
"plot_E_mk4.dat" using 1:2 with lines title "mk4"

set logscale y
set grid
plot "plot_E_esplicito.dat" using 1:4 with lines title "Eulero Esplicito", \
"plot_E_implicito.dat" using 1:4 with lines title "Eulero Implicito", \
"plot_Trapezio.dat" using 1:4 with lines title "Trapezi", \
"plot_E_mk4.dat" using 1:4 with lines title "mk4"