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
    for (double x = 2.f; x > 1e-20; x /= 2) {
        double result1 = log1p(x);
        double result2 = log(1 + x);
        printf("%10.5e %10.5e %10.5e\n", x, result1, result2);
    }
}