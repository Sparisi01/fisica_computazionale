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
#define lambda ((hc * hc * (M1 + M2)) / (2 * R * R * M1 * M2))
#define v ((V0) / (lambda))

// #define E_real -2.11990900
#define e 0.19040857383288309

#define A 1
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
        return A * sqrt(v - e) / (R * 2 * sqrt(PI));
    else if (x <= 1)
        return A * sin(sqrt(v - e) * x) / (x * R * 2 * sqrt(PI));
    else if (x > 1)
        return sin(sqrt(v - e)) / exp(-sqrt(e)) * exp(-x * sqrt(e)) / (x * R * 2 * sqrt(PI));
    return 0;
}

double psisquare(double x) {
    double p = psi(x);
    return p * p;
}

double funzione_integranda_numeratore(double x) {
    return pow(x, 2) * psisquare(x);
}

double funzione_integranda_numeratore2(double x) {
    return x * psisquare(x);
}

int main(int argc, char const *argv[]) {
    FILE *fileConvergenze;

    printf("%f\n", f(0.190409));
    fileConvergenze = fopen("convergenze.dat", "w");

    fprintf(fileConvergenze, "# Newton Secante Bisezione\n");

    double e_real = ricercaZeriBisezione(0.1, 3, 1000, f, 1.0e-11);
    double E_real = -e_real * lambda;
    for (size_t n_iterazioni = 1; n_iterazioni < 60; n_iterazioni++) {
        double e_bisezione = ricercaZeriBisezione(0.1, 3, n_iterazioni, f, 1.0e-11);
        double E_bisezione = -e_bisezione * lambda;

        double e_secante = ricercaZeriSecante(0.1, 3, n_iterazioni, f, 1.0e-11);
        double E_secante = -e_secante * lambda;

        double e_newton = ricercaZeriNewton(0.1, n_iterazioni, f, f_derivate, 1.0e-11);
        double E_newton = -e_newton * lambda;
        fprintf(fileConvergenze, "%d %10.15e %10.15e %10.15e\n", n_iterazioni, fabs(E_real - E_newton), fabs(E_real - E_secante), fabs(E_real - E_bisezione));
    }

    const char *commands = {
        "reset\n"
        "set datafile commentschars '#@&'\n"
        "set terminal png size 2048, 1536 font ', 36'\n"
        "set output 'plot.png'\n"
        "set key top right\n"
        "set grid \n"
        "set tics font ', 36'\n"
        "set title 'Convergenza metodi ricerca zeri' font ', 50'\n"
        "set xrange [0:40]\n"
        "set logscale y\n"
        "set xlabel 'Numero iterazioni'\n"
        "set ylabel '(E_{r} - E_{s}) [MeV]'\n"
        "plot 'convergenze.dat' using 1 : 4 with linespoint linewidth 3.5 pointsize 4 pointtype 6 title \"Bisezione\","
        "'' every ::0::7  using 1 : 3 with linespoint linewidth 3.5 pointsize 4 pointtype 6 title \"Secante\","
        "'' every ::0::5 using 1 : 2 with linespoint linewidth 3.5 pointsize 4 pointtype 6 title \"Newton\""};

    executeGNUPlotCommands(commands);

    const char *commands2 = {
        "reset\n"
        "set datafile commentschars '#@&'\n"
        "set terminal png size 2048, 1536 font ', 36'\n"
        "set output 'eqtrascendente.png'\n"
        "set key top right\n"
        "set grid \n"
        "set tics font ', 36'\n"
        "set title 'Stima soluzione equazione trascendente' font ', 50'\n"
        "set xrange [0:3]\n"
        "set xlabel 'e'\n"
        "set ylabel 'f(e)'\n"
        "v = 3.45804\n"
        "f(x)=1/tan(sqrt(v-x))+sqrt(x/(v-x))\n"
        "ax=0.190409\n"
        "ay=0.2\n"
        "plot f(x) with lines linewidth 3.5 linecolor 'blue' title 'Equazione (10)', 'soluzioneEnergia.dat' using 1 : 2 with points linewidth 3.5 pointsize 4 pointtype 6 linecolor 'blue' title 'Stima soluzione',"
        "1/tan(sqrt(v-x)) with lines linewidth 3.5 linecolor 'red' title 'Cotangente',"
        "-sqrt(x/(v-x)) with lines linewidth 3.5 linecolor 'green' title 'Radice' \n"};

    executeGNUPlotCommands(commands2);

    // Valor Medio

    double limsupint = 25;
    double norma = integraleSimpson(0, limsupint, 1000, psisquare);
    double media_posizione_quadra = integraleSimpson(0, limsupint, 1000, funzione_integranda_numeratore);
    double media_posizione = integraleSimpson(0, limsupint, 1000, funzione_integranda_numeratore2);

    double media_posizione_quadra_normalizzata = media_posizione_quadra / norma;
    double media_posizione_normalizzata = media_posizione / norma;

    double sigma = sqrt(media_posizione_quadra_normalizzata - media_posizione_normalizzata * media_posizione_normalizzata);

    printf("MEDIA POSIZIONE QUADRA: %f +- %f\n", media_posizione_normalizzata, sigma);

    FILE *plot_psi = fopen("valoriPsi.dat", "w");
    for (double i = 0; i < 6; i = i + 2.f / 1000.f) {
        fprintf(plot_psi, "%10.5e %10.5e\n", i, psisquare(i) / norma);
    }

    return 0;
}
