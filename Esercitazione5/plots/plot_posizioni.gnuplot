
set xrange [-10:10]
set yrange [-10:10]

plot "../data/posizione_init.dat" using 1:2 pt 6, \
"../data/posizione_end.dat" using 1:2 pt 6, \