#define PI 3.14159265359

#include <float.h>
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

    // ERRORE std distribuzione uniforme diviso per sqrt(N)

    // ESERCIZIO 3

    // 0.003383692573952728

    double sum_3 = 0;
    int M = 100000;

    for (size_t i = 0; i < M; i++)
    {
        double x = RND_Unif(0, 3);
        sum_3 += exp(-x * x / 2);
    }

    double I3a = sqrt(PI / 2) - sum_3 * 3 / M;
    printf("3.a: %.5E\n", I3a);

    // Si può usare direttamente la distribuzione gaussiana per fare la media
    // della funzione costante unitaria

    // # Math: <f> = \int_{-∞}^∞ θ(x-3)θ(x) \sqrt{\frac{1}{2\pi}} e^{-x^2/2}dx
    double sum_4 = 0;
    int K = 100000;

    for (size_t i = 0; i < K; i++)
    {
        double x = fabs(RND_Gauss(0, 1));
        sum_4 += (x <= 3 && x > 0) ? 1 : 0; // Distribuzione costante da 0 a 3
    }

    double I3b = sqrt(PI / 2) - sum_4 / K * sqrt(PI * 2) / 2;
    printf("3.b: %.5E\n", I3b);
    return 0;
}
