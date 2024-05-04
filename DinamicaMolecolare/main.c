#include "observables.h";
#include "system.h";
#include "verlet.h";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{

    // ALLOCAZIONE SISTEMA
    struct System s;
    s.n_particles = 10;
    s.t = 0;

    s.x = (double *)malloc(sizeof(double) * s.n_particles);
    s.y = (double *)malloc(sizeof(double) * s.n_particles);
    s.z = (double *)malloc(sizeof(double) * s.n_particles);

    if (!s.x || !s.y || s.z)
    {
        perror("Errore allocazione posizioni");
        exit(EXIT_FAILURE);
    }

    s.vx = (double *)malloc(sizeof(double) * s.n_particles);
    s.vy = (double *)malloc(sizeof(double) * s.n_particles);
    s.vz = (double *)malloc(sizeof(double) * s.n_particles);

    if (!s.vx || !s.vy || s.vz)
    {
        perror("Errore allocazione velocità");
        exit(EXIT_FAILURE);
    }

    s.masses = (double *)malloc(sizeof(double) * s.n_particles);

    if (!s.masses)
    {
        perror("Errore allocazione masse");
        exit(EXIT_FAILURE);
    }

    // VERLET
    const double t_in = 0;
    const double t_end = 1;
    const double dt = 1e-3;
    int n_dt = (t_end - t_in) / dt;

    for (size_t i = 0; i < n_dt; i++)
    {
    }

    return 0;
}
