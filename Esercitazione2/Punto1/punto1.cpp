#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S:\Coding\fisica_computazionale\utility\plot.h"
#include "S:\Coding\fisica_computazionale\utility\s_integrali.h"

#define alpha 2
#define beta 0.5
#define A 1
#define B 1

double f_a(double x) {
    double s = sin(alpha * x);
    return s * s;
}

double f_b(double x) {
    double s = sin(x);
    return x * x * s * s;
}

double f_c(double x) {
    return exp(-x * beta);
}

double f_d(double x) {
    return x * x * exp(-x * beta);
}

int main(int argc, char const* argv[]) {
    // Integrale a
    printf("(a) trapezi: %10.5e\n", integraleTrapezi(0, A, 100, f_a) - integraleTrapezi(0, A, 82, f_a));
    printf("(a) simpson: %10.5e\n", integraleSimpson(0, A, 100, f_a) - integraleSimpson(0, A, 10, f_a));

    // Integrale b
    printf("(b) trapezi: %10.5e\n", integraleTrapezi(0, A, 100, f_b) - integraleTrapezi(0, A, 61, f_b));
    printf("(b) simpson: %10.5e\n", integraleSimpson(0, A, 100, f_b) - integraleSimpson(0, A, 11, f_b));

    // Integrale c
    printf("(c) trapezi: %10.5e\n", integraleTrapezi(0, B, 100, f_c) - integraleTrapezi(0, B, 15, f_c));
    printf("(c) simpson: %10.5e\n", integraleSimpson(0, B, 100, f_c) - integraleSimpson(0, B, 2, f_c));

    // Integrale d
    printf("(d) trapezi: %10.5e\n", integraleTrapezi(0, B, 100, f_d) - integraleTrapezi(0, B, 41, f_d));
    printf("(d) simpson: %10.5e\n", integraleSimpson(0, B, 100, f_d) - integraleSimpson(0, B, 4, f_d));



    char const* gnuPlotCommands[] = {"plot \"data.tmp\""};

    return 0;
}
