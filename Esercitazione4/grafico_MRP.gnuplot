reset 
set multiplot layout 2,1 spacing 0.05
set monochrome
set title "MRP" font "Helvetica, 14"
set ylabel font ", 12"
set xlabel font ", 12"

set grid
set ytics nomirror
set xtics nomirror

set ylabel "Massa [M0]"

plot "./data/P_M_R_1.dat" u 5:4 w l notitle, \
"./data/P_M_R_2.dat" u 5:4 w l notitle, \
"./data/P_M_R_3.dat" u 5:4 w l notitle

set yrange [-0.5:3.5]
set xlabel "Raggio [R0]"
set notitle

set key right
set key box lt -1 lw 1
set key spacing 2 font "Helvetica, 12"

plot "./data/P_M_R_1.dat" u 3:2 w l t "stella 1", \
"./data/P_M_R_2.dat" u 3:2 w l t "stella 2",\
"./data/P_M_R_3.dat" u 3:2 w l t "stella 3"

