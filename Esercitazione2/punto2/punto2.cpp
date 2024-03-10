#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S:\Coding\fisica_computazionale\utility\plot.h"
#include "S:\Coding\fisica_computazionale\utility\s_integrali.h"
#include "S:\Coding\fisica_computazionale\utility\s_ricercaZeri.h"

#define V0 38.5  //[MeV]
// #define V0 27.4
#define R 1.93      //[fm]
#define hc 197.327  //[MeV]
#define M1 939.565  //[MeV]
#define M2 938.272  //[MeV]
#define lambda (hc * hc * (M1 + M2)) / (2 * R * R * M1 * M2)
#define v (V0) / (lambda)

#define E -2.11990900
#define e 0.190409

#define A 1
#define B 1
#define precisione 1e-7
#define PI 3.1415926535

using namespace std;

double f(double x) {
    return 1 / tan(sqrt(v - x)) + sqrt(x / (v - x));
}

double f_derivate(double x) {
    double sq = sqrt(v - x);
    double sin_sq = sin(sq);
    return 1 / (sin_sq * sin_sq * sq * 2) + v / (sqrt(x) * sq * (v - x) * 2);
}

// Funzioni integrale media posizione
double psi(double x) {
    if (x == 0)
        return A * sqrt(v - e) / (2 * sqrt(PI));
    else if (x <= R)
        return A * sin(sqrt(v - e) * x) / (x * 2 * sqrt(PI));
    else if (x > R)
        return B * exp(-x * sqrt(e)) / (x * 2 * sqrt(PI));

    return 0;
}

double psisquare(double x) {
    double p = psi(x);
    return p * p;
}

double funzione_integranda_numeratore(double x) {
    return pow(x, 4) * psisquare(x);
}

double funzione_integranda_numeratore2(double x) {
    return pow(x, 3) * psisquare(x);
}

int main(int argc, char const *argv[]) {
    FILE *fileConvergenze;

    fileConvergenze = fopen("convergenze.dat", "w");

    fprintf(fileConvergenze, "# Newton Secante Bisezione\n");

    for (size_t n_iterazioni = 1; n_iterazioni < 50; n_iterazioni++) {
        double e_bisezione = ricercaZeriBisezione(0.1, 3, n_iterazioni, f);
        double E_bisezione = -e_bisezione * lambda;

        double e_secante = ricercaZeriSecante(0.1, 3, n_iterazioni, f);
        double E_secante = -e_secante * lambda;

        double e_newton = ricercaZeriNewton(0.1, n_iterazioni, f, f_derivate);
        double E_newton = -e_newton * lambda;
        fprintf(fileConvergenze, "%d %10.5e %10.5e %10.5e\n", n_iterazioni, fabs(E - E_newton), fabs(E - E_secante), fabs(E - E_bisezione));
    }

    const char *commands = {
        "reset"
        "set datafile commentschars '#@&'\n"
        "set term png\n"
        "set output 'plot.png'\n"
        "set title 'Convergenza metodi ricerca zeri'\n"
        "set xrange [1:30]\n"
        "set logscale y\n"
        "set xlabel 'Numero iterazioni'\n"
        "set ylabel '(E_{r} - E_{s}) [MeV]'\n"
        "plot 'convergenze.dat' using 0 : 4 with linespoint title \"Bisezione\","
        "'' using 0 : 3 with linespoint title \"Secante\","
        "'' using 0 : 2 with linespoint title \"Newton\""};

    executeGNUPlotCommands(commands);

    // Valor Medio

    double limsupint = 30;
    double norma = integraleSimpson(0, limsupint, 1000, psisquare);
    double media_posizione_quadra = integraleSimpson(0, limsupint, 1000, funzione_integranda_numeratore);
    double media_posizione = integraleSimpson(0, limsupint, 1000, funzione_integranda_numeratore2);

    double media_posizione_quadra_normalizzata = media_posizione_quadra / norma;
    double media_posizione_normalizzata = media_posizione / norma;

    double sigma = sqrt(media_posizione_quadra_normalizzata - media_posizione_normalizzata * media_posizione_normalizzata);

    printf("MEDIA POSIZIONE QUADRA: %f +- %f\n", media_posizione_normalizzata, sigma);

    return 0;
}
