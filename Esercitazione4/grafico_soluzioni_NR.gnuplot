reset 
set multiplot layout 2,1 spacing 0.05
set title "Soluzione equazione di stabilita' stellare" font "Helvetica, 14"
set ylabel font ", 12"
set xlabel font ", 12"
set grid
set key right
set key box lt -1 lw 1
set key spacing 2 font "Helvetica, 12"
set grid
set ytics nomirror
set xtics nomirror
set format x ""
set lmargin 10

set ylabel "Pressione [P0]"
plot "./data/sol_stella_1_NR.dat" u 1:3 w l t "Fermionica NR", "./data/sol_stella_2_NR.dat" u 1:3 w l t "Fermionica UR", "./data/sol_stella_3_NR.dat" u 1:3 w l t "Nucleare"
unset title
unset format x
unset key
set xlabel "Raggio [R0]"
set ylabel "Massa [M0]"
plot "./data/sol_stella_1_NR.dat" u 1:2 w l notitle, "./data/sol_stella_2_NR.dat" u 1:2 w l notitle, "./data/sol_stella_3_NR.dat" u 1:2 w l notitle

unset multiplot



