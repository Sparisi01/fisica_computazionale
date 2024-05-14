set logscale x
set logscale y


plot "../data/monte_carlo_comparison.dat" using 1:2 with points t "Metodo 1", \
"" using 1:3 with points t "Metodo 2"
