reset 
set multiplot layout 1,2

set grid
set xlabel "Raggio"
set ylabel "Pressione"
plot "stella1.dat" using 1:3 with lines title "Pressione", "stella2.dat" using 1:3 with lines, "stella3.dat" using 1:3 with lines

set grid
set xlabel "Raggio"
set ylabel "Massa"
plot "stella1.dat" using 1:2 with lines title "Massa", "stella2.dat" using 1:2 with lines title "Massa", "stella3.dat" using 1:2 with lines title "Massa"

