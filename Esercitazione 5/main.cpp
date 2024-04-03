#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>

#define SIGMA 1
#define N 1
#define k 3
#define PI 3.1315926535897932384626433
using namespace std;

struct Vec3 {
    double x;
    double y;
    double z;
};

struct Particle {
    Vec3 *pos;
    Vec3 *vel;
    double massa;
};

struct System {
    double t;
    int N_particles;
    Particle *particles;
};

void verlet(System *system, double dt, Vec3 *(*F)(System *, double *), double *args) {
    Vec3 *oldForces = F(system, args);

    // Update position of all particles based on forces
    for (size_t i = 0; i < system->N_particles; i++) {
        Particle *p = &system->particles[i];
        p->pos->x += p->vel->x * dt + 0.5 / p->massa * oldForces[i].x * dt * dt;
        p->pos->y += p->vel->y * dt + 0.5 / p->massa * oldForces[i].y * dt * dt;
        p->pos->z += p->vel->z * dt + 0.5 / p->massa * oldForces[i].z * dt * dt;
    }
    // Update forces
    system->t += dt;
    Vec3 *newForces = F(system, args);
    // Update velocities
    for (size_t i = 0; i < system->N_particles; i++) {
        Particle *p = &system->particles[i];
        p->vel->x += 0.5 / p->massa * (oldForces[i].x + newForces[i].x) * dt;
        p->vel->y += 0.5 / p->massa * (oldForces[i].y + newForces[i].y) * dt;
        p->vel->z += 0.5 / p->massa * (oldForces[i].z + newForces[i].z) * dt;
    }
}

Vec3 *getForces(System *system, double *args) {
    struct Vec3 *forces = (Vec3 *)malloc(sizeof(struct Vec3) * system->N_particles);
    // Oscillatore Armonico
    for (size_t i = 0; i < system->N_particles; i++) {
        forces[i].x = -k * system->particles[i].pos->x;
        forces[i].y = 0;
        forces[i].z = 0;
    }

    return forces;
}

double *getRNDVelocity(double sigma) {
    double x[] = {rand() / (RAND_MAX + 1.),
                  rand() / (RAND_MAX + 1.)};

    double tmp;
    if (1 - x[0] == 1) {
        tmp = sigma * sqrt(-2 * log1p(1 - x[0]));
    } else {
        tmp = sigma * sqrt(-2 * log(1 - x[0]));
    }
    double y[] =
        {
            tmp * cos(2 * PI * x[1]),
            tmp * sin(2 * PI * x[1])};

    return y;
}

// Initialize velocities
void initVelocities(System *system, double T) {
    int n_velocities = system->N_particles * 3;
    bool is_heaven = (n_velocities % 2 == 0);

    if (is_heaven) {
        for (size_t i = 0; i < system->N_particles / 2; i += 2) {
            double *v1 = getRNDVelocity(T);
            double *v2 = getRNDVelocity(T);
            double *v3 = getRNDVelocity(T);
            system->particles[i].vel->x = v1[0];
            system->particles[i].vel->y = v1[1];
            system->particles[i].vel->z = v2[0];
            system->particles[i + 1].vel->x = v2[1];
            system->particles[i + 1].vel->y = v3[0];
            system->particles[i + 1].vel->z = v3[1];
        }
    } else {
    }
}

int main(int argc, char const *argv[]) {
    double t0 = 0;
    double dt = 1e-2;
    double tf = 10;
    double N_steps = (tf - t0) / dt;

    int seed = 8;

    System system;
    system.t = t0;
    system.N_particles = N;
    system.particles = (Particle *)malloc(sizeof(Particle) * N);

    // Inizialize masses
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[0].massa = 1;
    }

    // Inizialize velocities

    // Inizialize positions
    ofstream file("./data/evoluzione.dat");

    for (size_t i = 0; i < N_steps; i++) {
        verlet(&system, dt, getForces, NULL);

        file << system.t << " ";
        for (size_t i = 0; i < system.N_particles; i++) {
            file << system.particles[i].pos->x << " " << system.particles[i].pos->y << " " << system.particles[i].pos->z << " ";
            file << system.particles[i].vel->x << " " << system.particles[i].vel->y << " " << system.particles[i].vel->z << " ";
        }
        file << endl;
    }

    return 0;
}
