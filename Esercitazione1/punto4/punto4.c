#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    /**
     * Il peoblema è sempre lo stesso, quando sommo x a 1
     * se x non può essere traslato nei 52 bit del double
     * viene assorbito dall'1.
     * La funzione log1p sfrutta il calcolo polinomiale per ovviare
     * a questo problema difatti NON ho la somma.
     */

    FILE* file1;
    FILE* file2;
    FILE* file3;

    file1 = fopen("confronto_logaritmi.dat", "w");
    file2 = fopen("implementazione_stabile1.dat", "w");
    file3 = fopen("implementazione_stabile2.dat", "w");

    for (double x = 2; x > 1e-20; x /= 2) {
        // Confronto logaritmi
        double result1 = log1p(x);
        double result2 = log(1 + x);
        fprintf(file1, "%10.5e %10.5e %10.5e\n", x, result1, result2);

        // Stabile 1
        if (x + 1.f == 1.f) {
            fprintf(file2, "%10.5e\n", x);
        } else {
            double result = log(1.f + x);
            fprintf(file2, "%10.5e\n", result);
        }

        // Stabile 2
        if (x + 1.f == 1.f) {
            fprintf(file3, "%10.5e\n", x);
        } else {
            double result = x * (log(1.f + x) / ((1.f + x) - 1.f));
            fprintf(file3, "%10.5e\n", result);
        }
    }
}