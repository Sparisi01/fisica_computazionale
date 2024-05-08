#define PI 3.14159265359

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double RND_Unif(double a, double b)
{
    return (rand() / (RAND_MAX + 1.0)) * b - a;
}

double RND_Gauss(double mean, double sigma)
{
    double x = RND_Unif(0, 1);
    double y = RND_Unif(0, 1);

    return sigma * sqrt(-2 * log(1 - x)) * cos(2 * PI * y) + mean;
}

int main(int argc, char const *argv[])
{
    srand(1);

    // ESERCIZIO 1

    int N = 10000;

    double sum_1a = 0;
    double sum_1b = 0;
    double sum_1c = 0;
    double sum_1d = 0;

    for (size_t i = 0; i < N; i++)
    {
        double x = RND_Unif(0, PI);
        sum_1a += pow(sin(1 * x), 2);
        sum_1b += pow(x * sin(1 * x), 2);
        sum_1c += exp(-x);
        sum_1d += x * x * exp(-x);
    }

    double I1a = PI / N * sum_1a;
    double I1b = PI / N * sum_1b;
    double I1c = 1 - PI / N * sum_1c;
    double I1d = 2 - PI / N * sum_1d;

    printf("1.a: %.5E\n", I1a);
    printf("1.b: %.5E\n", I1b);
    printf("1.c: %.5E\n", I1c);
    printf("1.d: %.5E\n", I1d);

    return 0;
}
