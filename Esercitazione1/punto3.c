#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
    long double a = 1;
    long double c = 1;

    // I problemi nascono quando la differenza 4 * a * c + dell'ordine di errore_macchina * b*b
    for (long double b = 2; b < 1e12; b *= 2.) {
        long double x1 = (-b + sqrtl(b * b - 4 * a * c)) / (2 * a);
        long double x2 = (-b - sqrtl(b * b - 4 * a * c)) / (2 * a);

        long double x3 = 2 * c / (-b + sqrtl(b * b - 4 * a * c));
        long double x4 = 2 * c / (-b - sqrtl(b * b - 4 * a * c));
        printf("%10.5Le;%10.5Le\n", b * b, x1);
    }
}