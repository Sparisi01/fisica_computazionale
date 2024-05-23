#if !defined(MOTROPOLIS_H)
#define MOTROPOLIS_H

#include <stdlib.h>

double min(double a, double b)
{
    return a < b ? a : b;
}

double metropolis(double x, double (*w)(double y), double delta)
{
    double next_x = x + delta * rand() / (RAND_MAX + 1.); // Walker step
    double a = min(1, w(next_x) / w(x));                  // Calculate the a factor
    int accepted = rand() / (RAND_MAX + 1.) < a;

    return accepted ? next_x : x;
}

#endif // MOTROPOLIS_H
