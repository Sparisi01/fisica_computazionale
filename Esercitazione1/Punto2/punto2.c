#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S:\Coding\fisica_computazionale\utility\plot.h"

// gcc -o punto2.exe punto2.c -Wall -lm
// plot "dati.dat"

int main() {
    FILE* dati = fopen("dati_2.dat", "w");
    double e = exp(1.);
    for (double h = 2; h > 1e-15; h /= 2.) {
        double error = fabs((exp(1. + h) - e) / h - e);
        fprintf(dati, "%10.5e %10.5e\n", h, error);
    }

    const char* command =
        "reset\n"
        "set datafile commentschars '#@&'\n"
        "set terminal png size 2048, 1536 font 'Helvetica, 36'\n"
        "set output 'plot_2.png'\n"
        "set key top right\n"
        "set grid \n"
        "set title 'Errore in funzione del passo h' font ', 50'\n"
        "set xrange [1e-17:10]\n"
        "set yrange [1e-8:80]\n"
        "set logscale x\n"
        "set logscale y\n"
        "set format y '1x10^{%T}'\n"
        "set xlabel 'h'\n"
        "set ylabel 'Errore'\n"
        "e = exp(1)\n"
        "error = 2.2e-16\n"
        "f(x) = e*(error/x + x/2 + error + 5/4 * error*error*(1+1/x)+ 1/2*error*x + error*error*error*1/2/x)\n"
        "plot 'dati_2.dat' using 1 : 2 with points pointtype 6 pointsize 4 linewidth 3.5 linecolor 'blue' title 'Sperimentale',"
        "f(x) title 'Teorico' with lines linewidth 3.5 linecolor 'blue' ";

    executeGNUPlotCommands(command);
}
