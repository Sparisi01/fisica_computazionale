#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/** formatted output of ieee-754 representation of float */
void show_ieee754(float f) {
    union {
        float f;
        uint32_t u;
    } fu = {.f = f};
    int i = sizeof f * CHAR_BIT;

    printf("  ");
    while (i--)
        printf("%d ", (fu.u >> i) & 0x1);

    putchar('\n');
    printf(
        " |- - - - - - - - - - - - - - - - - - - - - - "
        "- - - - - - - - - -|\n");
    printf(
        " |s|      exp      |                  mantissa"
        "                   |\n\n");
}

// Metodo Stabile sia per b > 0  che eper b < 0
void radiciStabili(double a, double b, double c, FILE *file) {
    double x1, x2;

    if (b > 0) {
        x1 = (-b - sqrtl(b * b - 4 * a * c)) / (2 * a);
        x2 = 2 * c / (-b - sqrtl(b * b - 4 * a * c));
    } else {
        x1 = (-b + sqrtl(b * b - 4 * a * c)) / (2 * a);
        x2 = 2 * c / (-b + sqrtl(b * b - 4 * a * c));
    }

    fprintf(file, "%10.5e %10.5e %10.19e\n", b, x1, x2);
}

// Metodo instabile, fornisce x2 = 0
// per b^2 * epsilon = 4ac (x1 se b<0)
void radiciInstabili1(double a, double b, double c, FILE *file) {
    double x1 = (-b + sqrtl(-4 * a * c + b * b)) / (2 * a);
    double x2 = (-b - sqrtl(-4 * a * c + b * b)) / (2 * a);
    fprintf(file, "%10.5e %10.5e %10.5e %10.5e\n", b * b / (4 * a * c), x1, x2, b * b);
}

// Metodo instabile, fornisce x1 = inf per b^2 * epsilon = 4ac (x2 se b<0)
void radiciInstabili2(double a, double b, double c, FILE *file) {
    double x1 = 2 * c / (-b - sqrtl(b * b - 4 * a * c));
    double x2 = 2 * c / (-b + sqrtl(b * b - 4 * a * c));
    fprintf(file, "%10.5e %10.5e %10.5e\n", b * b / (4 * a * c), x1, x2);
}

int main() {
    double a = 1;
    double c = 1;

    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;
    // FILE *fptr4;

    fptr1 = fopen("stabile.dat", "w");
    fptr2 = fopen("instabile1.dat", "w");
    fptr3 = fopen("instabile2.dat", "w");
    // fptr4 = fopen("limite_instabilita.dat", "w");

    /**
     * 2^26 che è la radice di 2^52 il massimo numero esprimibile nella parte frazionaria
     * quando viene eseguita la somma la traslazione dell'sponenziale fa si che il numero
     * frazionario non possa più essere espresso in 52 bit e si ha quello che viene chiamato
     * l'effetto di assorbimento.
     */

    /* float n1 = 1;
    float n2 = pow(2, 24);
    show_ieee754(n1);
    show_ieee754(n2);
    show_ieee754(n1 + n2); */

    /*  radiciInstabili1(2, 4 * pow(2, 27), 2, fptr4);  // instabile

     radiciInstabili1(2, 4 * pow(2, 26), 2, fptr4);  // stabile */

    for (double b = 2; b < 1e15; b *= 2.) {
        radiciStabili(a, b, c, fptr1);
        // radiciStabili(a, -b, c, fptr1);
        radiciInstabili1(a, b, c, fptr2);
        // radiciInstabili1(a, -b, c, fptr2);
        radiciInstabili2(a, b, c, fptr3);
        // radiciInstabili2(a, -b, c, fptr3);
    }
}
