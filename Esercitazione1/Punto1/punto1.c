#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    //! L'errore si ha quando f è 2^n volte più piccolo di 2 dove n è il numero di bit dedicati alla frazione

    float f = 2;
    int n_div = 0;

    while (2.f + f != 2.f) {
        f = f / 2.f;
        n_div++;
    }

    printf("Errore float: %12.5e con %d divisioni\n", f, n_div);

    double d = 2;
    n_div = 0;
    while (2. + d != 2.) {
        d = d / 2.;
        n_div++;
    }
    printf("Errore double: %12.5e con %d divisioni\n", d, n_div);

    long double l = 2;
    n_div = 0;
    while (2.l + l != 2.l) {
        l = l / 2.l;
        n_div++;
    }
    printf("Errore long: %12.5Le con %d divisioni\n", l, n_div);
}