#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct vec3 {
    double tempo;
    double posizione;
    double velocita;
} vec3;

double H_oscillatore(vec3 q) {
    return 0.5 * (pow(q.posizione, 2) + pow(q.velocita, 2));
}

int euleroEsplicito(vec3* q, double t, int N, vec3 (*F)(vec3, double*), double* args, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    for (size_t i = 0; i < N; i++) {
        *q = (vec3){
            q->tempo + h,
            q->posizione + h * F(*q, args).velocita,
            q->velocita + h * F(*q, args).posizione,
        };
        fprintf(traiettoria, "%10.5e %10.5e %10.5e %10.5e\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 0;
}

int euleroEsplicitoStep(vec3* q, double h, vec3 (*F)(vec3, double*), double* args) {
    *q = (vec3){
        q->tempo + h,
        q->posizione + h * F(*q, args).velocita,
        q->velocita + h * F(*q, args).posizione,
    };
    return 1;
}

int euleroImplicitoOscillatore(vec3* q, double t, int N, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    if (h < 1e-8)
        return 0;
    for (size_t i = 0; i < N; i++) {
        *q = (vec3){
            q->tempo + h,
            (q->posizione + h * q->velocita) / (1 + h * h),
            (q->velocita - h * q->posizione) / (1 + h * h),
        };
        fprintf(traiettoria, "%10.5e %10.5e %10.5e %10.5e\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 1;
}

int euleroTrapeziOscillatore(vec3* q, double t, int N, FILE* traiettoria) {
    double h = (t - q->tempo) / N;
    if (h < 1e-8)
        return 0;
    for (size_t i = 0; i < N; i++) {
        *q = (vec3){
            q->tempo + h,
            q->posizione * (1 - h * h) + h * q->velocita,
            q->velocita - h * q->posizione,
        };
        fprintf(traiettoria, "%10.5e %10.5e %10.5e %10.5e\n", q->tempo, q->posizione, q->velocita, H_oscillatore(*q));
    }

    return 1;
}

vec3 F(vec3 q, double* parametri) {
    return (vec3){
        q.tempo,
        -q.posizione,
        q.velocita,
    };
}

int main(int argc, char const* argv[]) {
    FILE* plot_E_implicito = fopen("plot_E_implicito.dat", "w");
    FILE* plot_Trapezio = fopen("plot_Trapezio.dat", "w");
    FILE* plot_E_esplicito = fopen("plot_E_esplicito.dat", "w");

    vec3* vector1 = &(vec3){0, 1, 0};
    vec3* vector2 = &(vec3){0, 1, 0};
    vec3* vector3 = &(vec3){0, 1, 0};

    // euleroEsplicito(vector1, 250, 10000, F, NULL, plot_E_esplicito);
    double h = 250.f / 10000.f;
    for (size_t i = 0; i < 10000; i++) {
        euleroEsplicitoStep(vector1, h, F, NULL);
        fprintf(plot_E_esplicito, "%10.5e %10.5e %10.5e %10.5e\n", vector1->tempo, vector1->posizione, vector1->velocita, H_oscillatore(*vector1));
    }

    euleroImplicitoOscillatore(vector2, 250, 10000, plot_E_implicito);
    euleroTrapeziOscillatore(vector3, 250, 10000, plot_Trapezio);

    return 0;
}
