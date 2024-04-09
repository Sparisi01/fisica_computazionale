set grid
#splot "../data/starting_pos.dat" using 1:2:3
set xrange [-1:6]
set yrange [-1:6]
plot "../data/starting_pos.dat" using 2:3
