reset 

set ytics nomirror
set xtics nomirror
set xtics 10 
set xrange [0:100]
set multiplot layout 2,1 spacing 0.05 
set title "Massa e Raggio in funzione di P centrale" font "Helvetica, 18"

set lmargin 15
set grid
set ylabel "Massa [M0]" font "Helvetica, 12"

plot "./data/P_M_R_1.dat" u 5:4 w l notitle, \
"./data/P_M_R_2.dat" u 5:4 w l notitle, \
"./data/P_M_R_3.dat" u 5:4 w l notitle
set xrange [0:100]

set yrange [-0.5:3.5]
set xlabel "Raggio [R0]" font "Helvetica, 12"
set notitle

set key right box lt -1 lw 2 spacing 1.1 font "Helvetica, 14"
set bmargin 4

plot "./data/P_M_R_1.dat" u 3:2 w l t "Gamma = 5/3", \
"./data/P_M_R_2.dat" u 3:2 w l t "Gamma = 4/3",\
"./data/P_M_R_3.dat" u 3:2 w l t "Gamma = 2.54"

unset multiplot


