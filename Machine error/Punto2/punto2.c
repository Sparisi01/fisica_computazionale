#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S:\Coding\fisica_computazionale\utility\plot.h"

int main() {
    FILE* dati = fopen("dati_double.dat", "w");
    double e = exp(1.);
    for (double h = 2; h > 1e-15; h /= 2.) {
        double error = fabs((exp(1. + h) - e) / h - e);
        fprintf(dati, "%10.8e %10.8e\n", h, error);
    }

    FILE* dati2 = fopen("dati_float.dat", "w");
    float e_2 = expf(1.);
    for (float h_2 = 2; h_2 > 1e-8; h_2 /= 2.f) {
        float error_2 = fabsf((expf(1. + h_2) - e_2) / h_2 - e_2);
        fprintf(dati2, "%10.8e %10.8e\n", h_2, error_2);
    }
}
