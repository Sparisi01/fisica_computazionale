#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

    fprintf(file, "%10.5e %10.5e %10.5e\n", b, x1, x2);
}

// Metodo instabile, fornisce x2 = 0 per b^2 * epsilon = 4ac (x1 se b<0)
void radiciInstabili1(double a, double b, double c, FILE *file) {
    double x1 = (-b - sqrtl(b * b - 4 * a * c)) / (2 * a);
    double x2 = (-b + sqrtl(b * b - 4 * a * c)) / (2 * a);
    fprintf(file, "%10.5e %10.5e %10.5e\n", b, x1, x2);
}

// Metodo instabile, fornisce x1 = inf per b^2 * epsilon = 4ac (x2 se b<0)
void radiciInstabili2(double a, double b, double c, FILE *file) {
    double x1 = 2 * c / (-b + sqrtl(b * b - 4 * a * c));
    double x2 = 2 * c / (-b - sqrtl(b * b - 4 * a * c));
    fprintf(file, "%10.5e %10.5e %10.5e\n", b, x1, x2);
}

int main() {
    double a = 1;
    double c = 1;

    FILE *fptr1;
    FILE *fptr2;
    FILE *fptr3;

    fptr1 = fopen("stabile.dat", "w");
    fptr2 = fopen("instabile1.dat", "w");
    fptr3 = fopen("instabile2.dat", "w");

    for (double b = 2; b < 1e15; b *= 2.) {
        radiciStabili(a, b, c, fptr1);
        radiciStabili(a, -b, c, fptr1);
        radiciInstabili1(a, b, c, fptr2);
        radiciInstabili1(a, -b, c, fptr2);
        radiciInstabili2(a, b, c, fptr3);
        radiciInstabili2(a, -b, c, fptr3);
    }
}