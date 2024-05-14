#define PI 3.14159265359

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double RND_Unif(double min, double max)
{
    return min + rand() / (RAND_MAX + 1.0) * (max - min);
}

double RND_Gauss(double mean, double sigma)
{
    double x = RND_Unif(0, 1);
    double y = RND_Unif(0, 1);

    return mean + sigma * sqrt(-2 * log(1 - x)) * sin(2 * PI * y);
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

#define REAL_RESULT 0.003383692573952728

    FILE *monte_carlo_comparison = fopen("./data/monte_carlo_comparison.dat", "w");

    for (size_t i = 2; i < 7; i++)
    {
        double M = pow(8, i);

        double sum_3 = 0;
        double sum_4 = 0;

        for (size_t i = 0; i < M; i++)
        {
            double x = RND_Unif(0, 3);
            sum_3 += exp(-x * x / 2);

            // Si può usare direttamente la distribuzione gaussiana per fare la media
            // della funzione costante unitaria

            // # Math: <f> = 1/\sqrt{2\pi} \int_{-∞}^∞ θ(|x-3|)e^{-x^2/2} dx
            double y = RND_Gauss(0, 1);
            sum_4 += (fabs(y)) <= 3 ? 1 : 0; // Distribuzione costante da 0 a 3
        }

        double I3a = sqrt(PI / 2) - sum_3 * 3 / M;
        double I3b = sqrt(PI / 2) - sum_4 / M * sqrt(PI * 2) / 2;

        fprintf(monte_carlo_comparison, "%5.6E %5.6E %5.6E\n", M, fabs(I3a - REAL_RESULT), fabs(I3b - REAL_RESULT));
    }

    return 0;
}
