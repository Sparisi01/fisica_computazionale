#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double ricercaZeriSecante(double x0, double x1, int N, double f(double), double precision = 1e-7, bool verbose = false)
{
    if (f(x0) == 0) return x0;
    if (f(x1) == 0) return x1;

    double x2 = 0;
    int count = 1;
    do
    {
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        x0 = x1;
        x1 = x2;
        if (f(x2) == 0) return x2;
        if (verbose) printf("Iter %d: %10.5e\n", count, fabs(x0 - x1));
        count++;
    } while (fabs(x0 - x1) > precision && count < N);

    return x2;
}

double ricercaZeriNewton(double x0, int N, double f(double), double f_prime(double), double precision = 1e-7, bool verbose = false)
{
    if (f(x0) == 0) return x0;
    double x0_memory = 0;
    double x1 = 0;
    int count = 1;
    do
    {
        x1 = x0 - f(x0) / f_prime(x0);
        x0_memory = x0;
        x0 = x1;
        if (f(x1) == 0) return x1;
        if (verbose) printf("Iter %d: %10.5e\n", count, fabs(x0 - x1));
        count++;
    } while (fabs(x1 - x0_memory) > precision && count < N);

    return x0;
}

double ricercaZeriBisezione(double x0, double x1, int N, double f(double), double precision = 1e-7, bool verbose = false)
{
    if (f(x0) == 0) return x0;
    if (f(x1) == 0) return x1;
    if ((f(x0) < 0) - (f(x1) > 0) != 0) return NAN;

    double x2 = 0;
    int count = 1;

    do
    {
        x2 = x0 + (x1 - x0) / 2;
        if (f(x2) == 0) return x2;
        if (f(x2) > 0) x1 = x2;
        if (f(x2) < 0) x0 = x2;
        if (verbose) printf("Iter %d: %10.5e\n", count, fabs(x0 - x1));
        count++;
    } while (fabs(x0 - x1) > precision && count < N);

    return x0 + (x1 - x0) / 2;
}