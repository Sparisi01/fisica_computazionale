set terminal png size 2048, 1536 font ', 36' lw 4
set output "./out/grafico_soluzioni_R.png"
set multiplot layout 2,1 spacing 0.05
set title "Soluzione equazione di stabilità stellare relativistica" font ", 36"
set ylabel font ", 36"
set xlabel font ", 36"

set grid
set ytics nomirror
set xtics nomirror
set format x ""
set lmargin 10

set ylabel "Massa cumulativa [M_0]"
plot "./data/sol_stella_1_R.dat" u 1:2 w l notitle, "./data/sol_stella_2_R.dat" u 1:2 w l notitle, "./data/sol_stella_3_R.dat" u 1:2 w l notitle
unset title
unset format x
set key right
set key box lt -1 lw 1
set key spacing 2 font ", 36"
set xlabel "Raggio [R_0]"
set ylabel "Pressione [P_0]"
plot "./data/sol_stella_1_R.dat" u 1:3 w l t "Γ = 5/3 K = 0.05 RK", "./data/sol_stella_2_R.dat" u 1:3 w l t "Γ = 4/3 K = 0.1 RK", "./data/sol_stella_3_R.dat" u 1:3 w l t "Γ = 2.54 K = 0.01 RK"

unset multiplot



