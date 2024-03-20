reset 
set multiplot layout 1,2

set grid
set xlabel "Raggio"
set ylabel "Pressione"
plot "stella1.dat" using 1:3 with lines title "Stella 1", "stella2.dat" using 1:3 with lines title "Stella 2", "stella3.dat" using 1:3 with lines title "Stella 3"

set grid
set xlabel "Raggio"
set ylabel "Massa"
plot "stella1.dat" using 1:2 with lines title "Stella 1", "stella2.dat" using 1:2 with lines title "Stella 2", "stella3.dat" using 1:2 with lines title "Stella 3"

