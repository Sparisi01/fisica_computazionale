#include <math.h>
#include <stdio.h>

struct vecst {
    double t;
    double x;
    double v;
};

typedef struct vecst vecst;

int rk4(vecst* q, double h, vecst (*F)(vecst, double*), double* args) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    vecst tempF;

    k1 = h * F(*q, args).x;
    l1 = h * F(*q, args).v;

    tempF = F((vecst){q->t + h / 2, q->x + k1 / 2, q->v + l1 / 2}, args);
    k2 = h * tempF.x;
    l2 = h * tempF.v;

    tempF = F((vecst){q->t + h / 2, q->x + k2 / 2, q->v + l2 / 2}, args);
    k3 = h * tempF.x;
    l3 = h * tempF.v;

    tempF = F((vecst){q->t + h, q->x + k3, q->v + l3}, args);
    k4 = h * tempF.x;
    l4 = h * tempF.v;

    *q = (vecst){
        q->t + h,
        q->x + (k1 + 2 * (k2 + k3) + k4) / 6.,
        q->v + (l1 + 2 * (l2 + l3) + l4) / 6.,
    };

    return 1;
}

double ro(double P, double k, double gamma) {
    return pow(P / (k * (gamma - 1)), 1 / gamma);
}

vecst F(vecst q, double* parametri) {
    return (vecst){
        q.t,
        pow(q.t, 2) * ro(q.v, parametri[0], parametri[1]),
        -q.x * ro(q.v, parametri[0], parametri[1]) * pow(q.t, -2),
    };
}

int main(int argc, char const* argv[]) {
    FILE* stella1 = fopen("stella1.dat", "w");
    FILE* stella2 = fopen("stella2.dat", "w");
    FILE* stella3 = fopen("stella3.dat", "w");
    double T = 1;
    int N = 1000;
    double h = T / N;
    double Ps = 1;

    // k, gamma
    double parametri1[] = {0.05, 5. / 3.};
    double parametri2[] = {0.1, 4. / 3.};
    double parametri3[] = {0.01, 2.54};

    vecst* vector1 = &(vecst){0.001, 0, Ps};
    vecst* vector2 = &(vecst){0.001, 0, Ps};
    vecst* vector3 = &(vecst){0.001, 0, Ps};

    for (size_t j = 0; j < N; j++) {
        rk4(vector1, h, F, parametri1);
        rk4(vector2, h, F, parametri2);
        rk4(vector3, h, F, parametri3);
        fprintf(stella1, "%10.5E %10.5E %10.5E\n", vector1->t, vector1->x, vector1->v);
        fprintf(stella2, "%10.5E %10.5E %10.5E\n", vector2->t, vector2->x, vector2->v);
        fprintf(stella3, "%10.5E %10.5E %10.5E\n", vector3->t, vector3->x, vector3->v);
    }

    /* FILE* tante_stelle = fopen("tante_stelle.dat", "w");

    for (size_t i = 0; i < 20; i++)
    {
        double T = 1;

        int N = 1000;
        double h = T / N;
        double Ps = 1;
        vecst* vector = &(vecst){0.1, 0.01, Ps};
        for (size_t j = 0; j < N; j++) {
            rk4(vector, h, F, parametri1);
            fprintf(tante_stelle, "%10.5E %10.5E %10.5E\n", vector->t, vector->x, vector->v);
        }
}*/
    return 0;
}
