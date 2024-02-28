#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// gcc -o punto2.exe punto2.c -Wall -lm
// plot "dati.dat"

int main() {
    double e = exp(1.);
    for (double h = 1e-5; h > 1e-15; h /= 2.) {
        double error = fabs((exp(1. + h) - e) / h - e);
        printf("%10.5e %10.5e\n", h, error);
    }
}
