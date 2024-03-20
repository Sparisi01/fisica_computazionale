#include <math.h>
#include <stdio.h>
#include <stdlib.h>

struct vecst {
    double tempo;
    double posizione;
    double velocita;
};

typedef struct vecst vecst;

double H_oscillatore(vecst q) {
    return 0.5 * (pow(q.posizione, 2) + pow(q.velocita, 2));
}

int euleroEsplicito(vecst* q, double t, int N, vecst (*F)(vecst, double*), double* args, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    for (size_t i = 0; i < N; i++) {
        *q = (vecst){
            q->tempo + h,
            q->posizione + h * F(*q, args).velocita,
            q->velocita + h * F(*q, args).posizione,
        };
        fprintf(traiettoria, "%10.5E %10.5E %10.5E %10.5E\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 0;
}

int euleroEsplicitoStep(vecst* q, double h, vecst (*F)(vecst, double*), double* args) {
    *q = (vecst){
        q->tempo + h,
        q->posizione + h * F(*q, args).velocita,
        q->velocita + h * F(*q, args).posizione,
    };
    return 1;
}

int euleroImplicitoOscillatore(vecst* q, double t, int N, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    if (h < 1e-8)
        return 0;
    for (size_t i = 0; i < N; i++) {
        *q = (vecst){
            q->tempo + h,
            (q->posizione + h * q->velocita) / (1 + h * h),
            (q->velocita - h * q->posizione) / (1 + h * h),
        };
        fprintf(traiettoria, "%10.5E %10.5E %10.5E %10.5E\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 1;
}

int euleroTrapeziOscillatore(vecst* q, double t, int N, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    if (h < 1e-8)
        return 0;
    for (size_t i = 0; i < N; i++) {
        *q = (vecst){
            q->tempo + h,
            q->posizione * (1 - h * h) + h * q->velocita,
            q->velocita - h * q->posizione,
        };
        fprintf(traiettoria, "%10.5E %10.5E %10.5E %10.5E\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 1;
}

int rk4(vecst* q, double h, vecst (*F)(vecst, double*), double* args) {
    double k1, k2, k3, k4, l1, l2, l3, l4;
    vecst tempF;

    k1 = h * F(*q, args).posizione;
    l1 = h * F(*q, args).velocita;

    tempF = F((vecst){q->tempo + h / 2, q->posizione + k1 / 2, q->velocita + l1 / 2}, args);
    k2 = h * tempF.posizione;
    l2 = h * tempF.velocita;

    tempF = F((vecst){q->tempo + h / 2, q->posizione + k2 / 2, q->velocita + l2 / 2}, args);
    k3 = h * tempF.posizione;
    l3 = h * tempF.velocita;

    tempF = F((vecst){q->tempo + h, q->posizione + k3, q->velocita + l3}, args);
    k4 = h * tempF.posizione;
    l4 = h * tempF.velocita;

    *q = (vecst){
        q->tempo + h,
        q->posizione + (k1 + 2 * (k2 + k3) + k4) / 6,
        q->velocita + (l1 + 2 * (l2 + l3) + l4) / 6,
    };

    printf("%10.5E\n", q->velocita);

    return 1;
}

vecst F(vecst q, double* parametri) {
    return (vecst){
        q.tempo,
        q.velocita,
        -q.posizione,
    };
}

int main(int argc, char const* argv[]) {
    FILE* plot_E_implicito = fopen("plot_E_implicito.dat", "w");
    FILE* plot_Trapezio = fopen("plot_Trapezio.dat", "w");
    FILE* plot_E_esplicito = fopen("plot_E_esplicito.dat", "w");
    FILE* plot_E_mk4 = fopen("plot_E_mk4.dat", "w");

    vecst* vector1 = &(vecst){0, 1, 0};
    vecst* vector2 = &(vecst){0, 1, 0};
    vecst* vector3 = &(vecst){0, 1, 0};
    vecst* vector4 = &(vecst){0, 1, 0};

    double T = 25;
    int N = 1000;

    double h = T / N;
    for (size_t i = 0; i < N; i++) {
        euleroEsplicitoStep(vector1, h, F, NULL);
        fprintf(plot_E_esplicito, "%10.5E %10.5E %10.5E %10.5E\n", vector1->tempo, vector1->posizione, vector1->velocita, H_oscillatore(*vector1));
    }

    for (size_t j = 0; j < N; j++) {
        rk4(vector4, h, F, NULL);
        fprintf(plot_E_mk4, "%10.5E %10.5E %10.5E %10.5E\n", vector4->tempo, vector4->posizione, vector4->velocita, H_oscillatore(*vector4));
    }

    euleroImplicitoOscillatore(vector2, T, N, plot_E_implicito);
    euleroTrapeziOscillatore(vector3, T, N, plot_Trapezio);

    return 0;
}
