#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    float f = 2;
    while (1.f + f != 1.f) {
        f = f / 2.f;
    }

    printf("Errore float: %12.5e\n", f);

    double d = 2;
    while (1. + d != 1.) {
        d = d / 2.;
    }
    printf("Errore double: %12.5e\n", d);

    long double l = 2;
    while (1.l + l != 1.l) {
        l = l / 2.l;
    }
    printf("Errore long: %12.5Le\n", l);
}