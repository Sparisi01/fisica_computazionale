#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

struct Vec3 {
    double x, y, z;
};

struct Particle {
    Vec3 pos, vel;
    double mass;
};

struct System {
    double t;
    int N_particles;
    Particle *particles;
    Vec3 *forces;
};

#define PI 3.1315926535897932384626433  // To much accurate PI
#define SEED 9                          // Seed for random number generation
#define NEAREST_NEIGHBORS 1             // Choose how many boxes to keep in the lattice during interaction calculation
#define TEMP_IN 1.                      // Initial temperature
#define SIGMA_VELOCITIES sqrt(TEMP_IN)  //
#define N_PARTICLE 100                  //
#define SIGMA 1.                        // Distance scale
#define EPSILON 1.                      // Energy scale
#define DENSITY 1.2                     //
#define MASS 1.                         //
#define VOLUME (N_PARTICLE / DENSITY)   // Volume of a single cell (Volume in unit of sigma^3)
#define L cbrt(VOLUME)                  //
#define PRINT_THERMO 1                  //
#define PRINT_POSITIONS 0               //
#define PRINT_START_POSITIONS 0         //

void verlet(System &system, double dt, Vec3 *(*F)(const System &, double *), double *args) {
    Vec3 *oldForces = system.forces;
    // Update position of all particles based on forces
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].pos.x += system.particles[i].vel.x * dt + 0.5 / system.particles[i].mass * oldForces[i].x * dt * dt;
        system.particles[i].pos.y += system.particles[i].vel.y * dt + 0.5 / system.particles[i].mass * oldForces[i].y * dt * dt;
        system.particles[i].pos.z += system.particles[i].vel.z * dt + 0.5 / system.particles[i].mass * oldForces[i].z * dt * dt;
    }
    // Increase system time
    system.t += dt;
    // Update forces
    Vec3 *newForces = F(system, args);
    // Update velocities
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].vel.x += 0.5 / system.particles[i].mass * (oldForces[i].x + newForces[i].x) * dt;
        system.particles[i].vel.y += 0.5 / system.particles[i].mass * (oldForces[i].y + newForces[i].y) * dt;
        system.particles[i].vel.z += 0.5 / system.particles[i].mass * (oldForces[i].z + newForces[i].z) * dt;
    }

    // Update forces in system
    system.forces = newForces;
    free(oldForces);
}

double nearestDuplicateDistance(double x) {
    return x - L * rint(x / L);
}

double particlesSquareDistance(const Particle &p1, const Particle &p2) {
    double x = nearestDuplicateDistance(p1.pos.x - p2.pos.x);
    double y = nearestDuplicateDistance(p1.pos.y - p2.pos.y);
    double z = nearestDuplicateDistance(p1.pos.z - p2.pos.z);
    return pow(x, 2) + pow(y, 2) + pow(z, 2);
}

Vec3 particlesVersor(const Particle &p1, const Particle &p2) {
    double magnitude = sqrt(particlesSquareDistance(p1, p2));
    double x = nearestDuplicateDistance(p1.pos.x - p2.pos.x) / magnitude;
    double y = nearestDuplicateDistance(p1.pos.y - p2.pos.y) / magnitude;
    double z = nearestDuplicateDistance(p1.pos.z - p2.pos.z) / magnitude;
    return Vec3{.x = x, .y = y, .z = z};
}

Vec3 *getForcesOscillatore(const System &system, double *args) {
    struct Vec3 *forces = (Vec3 *)malloc(sizeof(struct Vec3) * system.N_particles);
    double k1 = args[0];
    double k2 = args[1];
    double k3 = args[2];
    for (size_t i = 0; i < system.N_particles; i++) {
        forces[i].x = -k1 * system.particles[i].pos.x;
        forces[i].y = -k2 * system.particles[i].pos.y;
        forces[i].z = -k3 * system.particles[i].pos.z;
    }
    return forces;
}

Vec3 *getForcesLennarJones(const System &system, double *args) {
    struct Vec3 *forces = (Vec3 *)calloc(system.N_particles, sizeof(struct Vec3));
    for (size_t i = 0; i < system.N_particles - 1; i++) {
        for (size_t j = i + 1; j < system.N_particles; j++) {
            // --------------------------------
            double r_square = particlesSquareDistance(system.particles[i], system.particles[j]);
            double r = sqrt(r_square);
            Vec3 r_dir_ij = particlesVersor(system.particles[i], system.particles[j]);
            double F_magnitude = 24 * EPSILON * SIGMA / r_square * (2 * pow(SIGMA / r, 11) - pow(SIGMA / r, 5));
            Vec3 force_ij = {r_dir_ij.x * F_magnitude, r_dir_ij.y * F_magnitude, r_dir_ij.z * F_magnitude};

            // --------------------------------
            forces[i].x += force_ij.x;
            forces[i].y += force_ij.y;
            forces[i].z += force_ij.z;
            // Newton third principle
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
void initVelocities(System &system, double sigma) {
    for (size_t i = 0; i < system.N_particles; i += 1) {
        system.particles[i].vel.x = getRNDVelocity(sigma);
        system.particles[i].vel.y = getRNDVelocity(sigma);
        system.particles[i].vel.z = getRNDVelocity(sigma);
    }
}

// Initialize positions in a cubic lattice
void initPositions(System &system, int N_particle) {
    double N_seg_particle = ceil(cbrt(N_particle));

    // Correzione per lasciare metà delle superfici vuote sul cubo in
    // modo da evitare la sovrapposizione tra le particelle e provocare
    // una distanza reciproca nulla
    double jump = L / (N_seg_particle + 1);

    int q = 0;
    for (size_t i = 0; i < N_seg_particle; i++) {
        for (size_t j = 0; j < N_seg_particle; j++) {
            for (size_t k = 0; k < N_seg_particle; k++) {
                if (q == N_particle) {
                    return;
                } else {
                    system.particles[q].pos = Vec3{.x = i * jump, .y = j * jump, .z = k * jump};
                    q++;
                }
            }
        }
    }
}

double mean_array(const double *array, int len) {
    double sum = 0;
    for (size_t i = 0; i < len; i++) {
        sum += array[i];
    }

    return sum / len;
}

double var_array(const double *array, int len) {
    double mean = mean_array(array, len);
    double sum = 0;
    for (size_t i = 0; i < len; i++) {
        sum += pow(array[i] - mean, 2);
    }

    return sum / len;
}

double TQ_Temperature(const System &system) {
    double sum = 0;
    for (size_t i = 0; i < system.N_particles; i++) {
        sum += system.particles[i].mass * pow(system.particles[i].vel.x, 2);
        sum += system.particles[i].mass * pow(system.particles[i].vel.y, 2);
        sum += system.particles[i].mass * pow(system.particles[i].vel.z, 2);
    }
    return sum / (3 * system.N_particles);
}

double TQ_Pressure(const System &system) {
    return 0;
}

void printCMVelocity(const System &system, FILE *file = stdout) {
    // Verifica pos centro di massa
    double CM_x = 0;
    double CM_y = 0;
    double CM_z = 0;

    for (size_t i = 0; i < system.N_particles; i++) {
        CM_x += system.particles[i].vel.x;
        CM_y += system.particles[i].vel.y;
        CM_z += system.particles[i].vel.z;
    }

    fprintf(stdout, "%10.15E %10.15E %10.15E\n", CM_x / system.N_particles, CM_y / system.N_particles, CM_z / system.N_particles);
}

int main(int argc, char const *argv[]) {
    double t0 = 0.;    // Time zero simulation
    double tf = 1.;    // Time final simulation
    double dt = 1e-4;  // Time interval
    double N_time_steps = (tf - t0) / dt;
    srand(SEED);
    // FILES
    FILE *starting_postions_file = fopen("./data/starting_pos.dat", "w");
    FILE *particle_data_file = fopen("./data/particles_data.dat", "w");
    FILE *termo_data_file = fopen("./data/termo_data.dat", "w");

    // SYSTEM INITIALIZATION
    System system;
    system.t = t0;
    system.N_particles = N_PARTICLE;
    system.particles = (Particle *)malloc(sizeof(Particle) * system.N_particles);

    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].mass = MASS;
    }

    initVelocities(system, SIGMA_VELOCITIES);

    initPositions(system, system.N_particles);

    system.forces = getForcesLennarJones(system, NULL);

    if (PRINT_START_POSITIONS) {
        for (size_t i = 0; i < system.N_particles; i++) {
            fprintf(starting_postions_file, "%f %f %f %f %f %f \n", system.particles[i].pos.x, system.particles[i].pos.y, system.particles[i].pos.z, system.particles[i].vel.x, system.particles[i].vel.y, system.particles[i].vel.z);
        }
    }

    // Initialize thermodinamic variables arrays
    double *temperature_array = (double *)malloc(sizeof(double) * N_time_steps);
    double *pressure_array = (double *)malloc(sizeof(double) * N_time_steps);

    // System evolution
    for (size_t i = 0; i < N_time_steps; i++) {
        // printCMVelocity(&system);

        // Calculate thermodinamics variables
        temperature_array[i] = TQ_Temperature(system);
        pressure_array[i] = TQ_Pressure(system);

        if (PRINT_POSITIONS) {
            // Print particles info
            fprintf(particle_data_file, "%f ", system.t);
            for (size_t i = 0; i < system.N_particles; i++) {
                fprintf(particle_data_file, "%f %f %f %f %f %f ", system.particles[i].pos.x, system.particles[i].pos.y, system.particles[i].pos.z, system.particles[i].vel.x, system.particles[i].vel.y, system.particles[i].vel.z);
            }
            fprintf(particle_data_file, "\n");
        }

        verlet(system, dt, getForcesLennarJones, NULL);
    }

    if (PRINT_THERMO) {
        for (size_t i = 0; i < N_time_steps; i++) {
            fprintf(termo_data_file, "%f %f %f\n", t0 + i * dt, temperature_array[i], pressure_array[i]);
        }

        printf("Temperatura: %f +- %f\n", mean_array(temperature_array, N_time_steps), sqrt(var_array(temperature_array, N_time_steps)));
        printf("Pressione: %f +- %f\n", mean_array(pressure_array, N_time_steps), sqrt(var_array(pressure_array, N_time_steps)));
    }

    free(temperature_array);
    free(pressure_array);
    return 0;
}
