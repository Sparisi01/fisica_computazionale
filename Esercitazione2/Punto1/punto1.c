#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double f(double x);
double metodoTrapezi(double a, double b, int N);
double metodoSimpson(double a, double b, int N);

int main(int argc, char const *argv[]) {
    printf("%10.5e\n", metodoTrapezi(0, 1, 1000));
    printf("%10.5e\n", metodoSimpson(0, 1, 1000));
    return 0;
}

double metodoTrapezi(double a, double b, int N) {
    double h = (b - a) / N;
    double sum = h / 2 * (f(a) + f(b));

    for (size_t i = 1; i < N; i++) {
        sum += h * f(a + i * h);
    }

    return sum;
}

double metodoSimpson(double a, double b, int N) {
    if (N % 2 == 1) N++;
    double h = (b - a) / N;
    double sum = f(a) + f(b);

    for (size_t i = 1; i <= N / 2 - 1; i++) {
        sum += 2 * f(a + 2 * i * h) + 4 * f(a + ((2 * i) - 1) * h);
    }

    return (sum + 4 * f(a + (N - 1) * h)) * h / 3;
}

double f(double x) {
    return exp(x);
}