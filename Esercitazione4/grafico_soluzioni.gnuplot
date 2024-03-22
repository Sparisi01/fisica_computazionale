reset 

set monochrome
set multiplot layout 2,1 spacing 0.05
set title "Soluzione equazione di stabilità stellare" font "Helvetica, 14"
set ylabel font "Helvetica, 12"
set xlabel font "Helvetica, 12"
set key right
set key box lt -1 lw 1
set key spacing 2 font "Helvetica, 12"
set grid
set ytics nomirror
set xtics nomirror
set format x ""

set ylabel "Pressione [P0]"
plot "stella1.dat" u 1:3 w l t "Fermionica NR", "stella2.dat" u 1:3 w l t "Fermionica UR", "stella3.dat" u 1:3 w l t "Nucleare"
unset title
unset format x
unset key
set xlabel "Raggio [R0]"
set ylabel "Massa [M0]"
plot "stella1.dat" u 1:2 w l notitle, "stella2.dat" u 1:2 w l notitle, "stella3.dat" u 1:2 w l notitle

