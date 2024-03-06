#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "S:\Coding\fisica_computazionale\algoritmi\s_integrali.h"

double f(double x) {
    return exp(x);
}

int main(int argc, char const *argv[]) {
    printf("%10.5e\n", integraleTrapezi(0, 1, 1000, f));
    printf("%10.5e\n", integraleSimpson(0, 1, 1000, f));
    return 0;
}
