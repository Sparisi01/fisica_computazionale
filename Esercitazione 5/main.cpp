#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>

#define PI 3.1315926535897932384626433
using namespace std;

struct Vec3 {
    double x;
    double y;
    double z;
};

struct Particle {
    Vec3 pos;
    Vec3 vel;
    double mass;
};

struct System {
    double t;
    int N_particles;
    Particle *particles;
    Vec3 *forces;
};

void verlet(System *system, double dt, Vec3 *(*F)(System *, double *), double *args) {
    Vec3 *oldForces = system->forces;
    // Update position of all particles based on forces
    for (size_t i = 0; i < system->N_particles; i++) {
        system->particles[i].pos.x += system->particles[i].vel.x * dt + 0.5 / system->particles[i].mass * oldForces[i].x * dt * dt;
        system->particles[i].pos.y += system->particles[i].vel.y * dt + 0.5 / system->particles[i].mass * oldForces[i].y * dt * dt;
        system->particles[i].pos.z += system->particles[i].vel.z * dt + 0.5 / system->particles[i].mass * oldForces[i].z * dt * dt;
    }
    // Increase system time
    system->t += dt;
    // Update forces
    Vec3 *newForces = F(system, args);
    // Update velocities
    for (size_t i = 0; i < system->N_particles; i++) {
        system->particles[i].vel.x += 0.5 / system->particles[i].mass * (oldForces[i].x + newForces[i].x) * dt;
        system->particles[i].vel.y += 0.5 / system->particles[i].mass * (oldForces[i].y + newForces[i].y) * dt;
        system->particles[i].vel.z += 0.5 / system->particles[i].mass * (oldForces[i].z + newForces[i].z) * dt;
    }

    // Update forces in system
    system->forces = newForces;
}

Vec3 *getForces(System *system, double *args) {
    struct Vec3 *forces = (Vec3 *)malloc(sizeof(struct Vec3) * system->N_particles);

    /* // Oscillatore Armonico
    for (size_t i = 0; i < system->N_particles; i++) {
        forces[i].x = -k * system->particles[i].pos.x;
        forces[i].y = 0;
        forces[i].z = 0;
    } */

    for (size_t i = 0; i < system->N_particles - 1; i++) {
        for (size_t j = i + 1; j < system->N_particles; j++) {
            Vec3 force_ij = NULL;
            forces[i].x += force_ij.x;
            forces[i].y += force_ij.y;
            forces[i].z += force_ij.z;
            // 3 principle
            forces[j].x -= force_ij.x;
            forces[j].y -= force_ij.y;
            forces[j].z -= force_ij.z;
        }
    }

    return forces;
}

double getRNDVelocity(double sigma) {
    double x = rand() / (RAND_MAX + 1.);

    double tmp;
    if (1 - x == 1) {
        tmp = sigma * sqrt(-2 * log1p(x));
    } else {
        tmp = sigma * sqrt(-2 * log(1 - x));
    }
    return tmp * cos(2 * PI * x);
}

// Initialize velocities
void initVelocities(System *system, double sigma) {
    for (size_t i = 0; i < system->N_particles; i += 1) {
        system->particles[i].vel.x = getRNDVelocity(sigma);
        system->particles[i].vel.y = getRNDVelocity(sigma);
        system->particles[i].vel.z = getRNDVelocity(sigma);
    }
}

// Initialize positions
void initPositions(System *system, double sigma) {
    for (size_t i = 0; i < system->N_particles; i += 1) {
        system->particles[i].pos.x = getRNDVelocity(sigma);
        system->particles[i].pos.y = getRNDVelocity(sigma);
        system->particles[i].pos.z = getRNDVelocity(sigma);
    }
}

double TQ_Temperature(System *system) {
    double sum = 0;
    for (size_t i = 0; i < system->N_particles; i++) {
        sum += system->particles[i].mass * pow(system->particles[i].vel.x, 2);
        sum += system->particles[i].mass * pow(system->particles[i].vel.y, 2);
        sum += system->particles[i].mass * pow(system->particles[i].vel.z, 2);
    }
    return sum / (3 * system->N_particles);
}

double TQ_Pressure(System *system) {
}

int main(int argc, char const *argv[]) {
    double const t0 = 0;          // Time zero simulation
    double const tf = 10;         // Time final simulation
    double const dt = 1e-2;       // Time interval
    double const Temp_in = 1;     // Initial temperature
    int const seed = 8;           // Seed for random number generation
    double const sigma_vel = 1;   // Sigma velocity distribution
    double const N_particle = 1;  // Number of particles in simulation

    double const N_time_steps = (tf - t0) / dt;
    srand(seed);

    System system;
    system.t = t0;
    system.N_particles = N_particle;
    system.particles = (Particle *)malloc(sizeof(Particle) * system.N_particles);

    // Inizialize masses
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].mass = 1;
        system.particles[i].vel = {.x = 1, .y = 1, .z = 1};
        system.particles[i].pos = {.x = 1, .y = 1, .z = 1};
    }

    // Init velocities
    // initVelocities(&system, sigma_vel);

    // Init positions
    // initPositions(&system, 0);

    // Init Forces
    system.forces = getForces(&system, NULL);

    ofstream file("./data/evoluzione.dat");

    for (size_t i = 0; i < N_time_steps; i++) {
        verlet(&system, dt, getForces, NULL);

        file << system.t << " ";
        for (size_t i = 0; i < system.N_particles; i++) {
            file << system.particles[i].pos.x << " " << system.particles[i].pos.y << " " << system.particles[i].pos.z << " ";
            file << system.particles[i].vel.x << " " << system.particles[i].vel.y << " " << system.particles[i].vel.z << " ";
        }
        file << endl;
    }

    return 0;
}
