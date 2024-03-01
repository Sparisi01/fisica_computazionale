epsi =  1.11022e-16
set logscale x
set logscale y
set xrange [1e-15:3]
set grid
plot exp(1)*(epsi/x + x/2), "dati_2.dat" title "Errore derivata Double"
