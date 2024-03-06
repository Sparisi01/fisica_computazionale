#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define v 3.4584
#define R 1.93      //[fm]
#define hc 197.327  //[MeV]
#define M1 939.565  //[MeV]
#define M2 938.272  //[MeV]

long double f(long double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}

long double secante(long double x0, long double x1, int N) {
    long double x2 = 0;
    int count = 1;
    do {
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
        printf("Iter %d: %10.5e\n", count, fabs(x0 - x1));
        count++;
    } while (fabs(x0 - x1) > 1e-7 && count < N);

    return x2;
}

long double newton(long double x0, long double x1, int N) {
    long double x2 = 0;
    int count = 1;
    do {
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
        printf("Iter %d: %10.5e\n", count, fabs(x0 - x1));
        count++;
    } while (fabs(x0 - x1) > 1e-12 && count < N);

    return x2;
}

int main(int argc, char const *argv[]) {
    long double e = secante(2, 3, 100);
    long double lambda = hc * hc * (M1 + M2) / (2 * R * R * M1 * M2);
    long double E = -e * lambda;
    printf("e = %10.5Le\nE = %10.5Le\n", e, E);
    return 0;
}
