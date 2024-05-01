#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <iostream>

using namespace std;
using namespace std::chrono;
struct Vec3 {
    double x, y, z;

    void operator+=(const Vec3 &a) {
        this->x += a.x;
        this->y += a.y;
        this->z += a.z;
    }

    void operator-=(const Vec3 &a) {
        this->x -= a.x;
        this->y -= a.y;
        this->z -= a.z;
    }

    // Dot product
    double operator*(const Vec3 &a) const {
        return a.x * x + a.y * y + a.z * z;
    }
    // Const product
    Vec3 operator*(const double a) const {
        return Vec3{a * x, a * y, a * z};
    }

    Vec3 operator/(const double a) const {
        return Vec3{x / a, y / a, z / a};
    }
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
    double kinetic_en;
    double pot_en;
    double *forcesWork;
    double densita;
    double volume;
    double L;
};

struct InitialCondition {
    double densita;
    double temperatura;
    int stampa;
    const char *file_name_g;
    const char *file_name_thermo;
};
// TEMP - RHO - M - N_CELLS - TIME_CUT - TEq
// 2 - 1.2 - 4 - 4 - 0.4 - 1.062724

#define PI 3.14159265359
#define SEED 9

// UNITÀ DI MISURA
#define SIGMA 1.
#define EPSILON 1.
#define MASS 1.
// STRUTTURA RETICOLO
#define M 4  // M=1 CC, M=2 BCC, M=4 FCC
#define N_CELLS 5
#define FREQ 0
// OUTPUT
#define PRINT_THERMO 1

void verletPropagator(System &system, double dt, Vec3 *(*F)(System &, double *), double *args) {
    system.kinetic_en = 0;
    system.pot_en = 0;
    Vec3 *oldForces = system.forces;
    // Aggiorna posizioni di tutte le particelle
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].pos.x += system.particles[i].vel.x * dt + 0.5 / system.particles[i].mass * oldForces[i].x * dt * dt;
        system.particles[i].pos.y += system.particles[i].vel.y * dt + 0.5 / system.particles[i].mass * oldForces[i].y * dt * dt;
        system.particles[i].pos.z += system.particles[i].vel.z * dt + 0.5 / system.particles[i].mass * oldForces[i].z * dt * dt;
    }
    // Incrementa il tempo del sistema
    system.t += dt;
    // Calcola forze agenti sulle particelle
    auto tin = high_resolution_clock::now();
    Vec3 *newForces = F(system, args);
    auto tend = high_resolution_clock::now();
    // printf("%d ms \n", duration_cast<milliseconds>(tend - tin));
    //  Aggiorna velocità
    for (size_t i = 0; i < system.N_particles; i++) {
        system.kinetic_en += (system.particles[i].vel * system.particles[i].vel) / (2 * system.particles[i].mass);

        system.particles[i].vel.x += 0.5 / system.particles[i].mass * (oldForces[i].x + newForces[i].x) * dt;
        system.particles[i].vel.y += 0.5 / system.particles[i].mass * (oldForces[i].y + newForces[i].y) * dt;
        system.particles[i].vel.z += 0.5 / system.particles[i].mass * (oldForces[i].z + newForces[i].z) * dt;
    }

    system.forces = newForces;
    free(oldForces);
    // oldForces = NULL;
}

Vec3 *getForcesOscillatore(System &system, double *args) {
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

double lennarJonesPotential(const System &system, double r) {
    if (r >= system.L / 2) {
        return 0;
    } else {
        double VSHIFT = (4 * EPSILON * (pow(SIGMA / (system.L / 2), 12) - pow(SIGMA / (system.L / 2), 6)));
        return 4 * EPSILON * (pow(SIGMA / r, 12) - pow(SIGMA / r, 6)) - VSHIFT;
    }
}

Vec3 *getForcesLennarJones(System &system, double *args) {
    struct Vec3 *forces = (Vec3 *)calloc(system.N_particles, sizeof(struct Vec3));
    // Reset lavori forze
    for (size_t i = 0; i < system.N_particles; i++)
        system.forcesWork[i] = 0;

    // Cicli for per ogni particella
    for (size_t i = 0; i < system.N_particles - 1; i++) {
        for (size_t j = i + 1; j < system.N_particles; j++) {
            Vec3 &p1 = system.particles[i].pos;
            Vec3 &p2 = system.particles[j].pos;
            double x = (p1.x - p2.x) - system.L * rint((p1.x - p2.x) / system.L);
            double y = (p1.y - p2.y) - system.L * rint((p1.y - p2.y) / system.L);
            double z = (p1.z - p2.z) - system.L * rint((p1.z - p2.z) / system.L);
            double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
            Vec3 r_dir_ij = Vec3{.x = x, .y = y, .z = z} / r;
            // double F_magnitude = 24 * EPSILON * SIGMA / r_square * (2 * pow(SIGMA / r, 11) - pow(SIGMA / r, 5));
            double h = 1e-8;
            system.pot_en += lennarJonesPotential(system, r);
            double F_magnitude = -(lennarJonesPotential(system, r + h) - lennarJonesPotential(system, r - h)) / (2 * h);
            Vec3 force_ij = r_dir_ij * (F_magnitude);
            forces[i] += force_ij;
            forces[j] -= force_ij;

            system.forcesWork[i] += force_ij * (r_dir_ij * r);
            system.forcesWork[j] += force_ij * (r_dir_ij * r);
        }
    }

    return forces;
}

double getRNDVelocity(double SIGMA_VELOCITIES) {
    double x = rand() / (RAND_MAX + 1.);

    double tmp;
    if (1 - x == 1) {
        tmp = SIGMA_VELOCITIES * sqrt(-2 * log1p(x));
    } else {
        tmp = SIGMA_VELOCITIES * sqrt(-2 * log(1 - x));
    }
    return tmp * cos(2 * PI * x);
}

// Initialize velocities
void initVelocities(System &system, double SIGMA_VELOCITIES) {
    for (size_t i = 0; i < system.N_particles; i += 1) {
        system.particles[i].vel.x = getRNDVelocity(SIGMA_VELOCITIES);
        system.particles[i].vel.y = getRNDVelocity(SIGMA_VELOCITIES);
        system.particles[i].vel.z = getRNDVelocity(SIGMA_VELOCITIES);
    }
}

// Initialize positions in a cubic lattice
void initPositions(System &system) {
    Vec3 lattice_positions[M];
    switch (M) {
        case 1:  // CC
            lattice_positions[0] = {.x = 0, .y = 0, .z = 0};
            break;
        case 2:  // BCC
            lattice_positions[0] = {.x = 0, .y = 0, .z = 0};
            lattice_positions[1] = {.x = 0.5, .y = 0.5, .z = 0.5};
            break;
        case 4:  // FCC
            lattice_positions[0] = {.x = 0, .y = 0, .z = 0};
            lattice_positions[1] = {.x = 0, .y = 0.5, .z = 0.5};
            lattice_positions[2] = {.x = 0.5, .y = 0.5, .z = 0};
            lattice_positions[3] = {.x = 0.5, .y = 0, .z = 0.5};
            break;
        default:
            break;
    }

    double jump = system.L / N_CELLS;
    int q = 0;
    for (size_t i = 0; i < N_CELLS; i++) {
        for (size_t j = 0; j < N_CELLS; j++) {
            for (size_t k = 0; k < N_CELLS; k++) {
                for (size_t w = 0; w < M; w++) {
                    system.particles[q].pos = Vec3{.x = (i + lattice_positions[w].x), .y = (j + lattice_positions[w].y), .z = (k + lattice_positions[w].z)} * jump;
                    q++;
                }
            }
        }
    }
}

double mean_array(const double *array, int len, int first) {
    double sum = 0;
    for (size_t i = first; i < len; i++) {
        sum += array[i];
    }

    return sum / (len - first);
}

double var_array(const double *array, int len, int first) {
    double mean = mean_array(array, len, first);
    double sum = 0;
    for (size_t i = first; i < len; i++) {
        sum += pow(array[i] - mean, 2);
    }

    return sum / (len - first);
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
    double W = 0;
    for (size_t i = 0; i < system.N_particles; i++) {
        W += system.forcesWork[i];
    }
    W /= 2 * system.N_particles;
    double T = TQ_Temperature(system);

    // return (T * system.densita + W / 3);  // Pressione
    return (1 + W / (3 * T));  // Compressibiità
}

void distribuzioneRadiale(const System &system, FILE *file) {
    double max_radius = system.L / 2 * sqrt(3);
    double N_intervals = 150;
    double radius_interval = max_radius / N_intervals;
    int *counting_array = (int *)calloc(sizeof(int), N_intervals);
    // -----------------------------------
    for (size_t j = 0; j < system.N_particles; j++) {
        Vec3 &p1 = system.particles[j].pos;
        for (size_t i = 0; i < system.N_particles; i++) {
            if (i == j) continue;

            double x = (p1.x - system.particles[i].pos.x) - system.L * rint((p1.x - system.particles[i].pos.x) / system.L);
            double y = (p1.y - system.particles[i].pos.y) - system.L * rint((p1.y - system.particles[i].pos.y) / system.L);
            double z = (p1.z - system.particles[i].pos.z) - system.L * rint((p1.z - system.particles[i].pos.z) / system.L);
            double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
            int n_radius_jump = (int)round(r / radius_interval);
            counting_array[n_radius_jump]++;
        }
    }
    for (size_t i = 0; i < N_intervals; i++) {
        fprintf(file, "%f %f\n", i * radius_interval * 2 / system.L, (double)(counting_array[i]) / (system.densita * system.N_particles * (4 * PI / 3) * (pow((i + 1) * radius_interval, 3) - pow(i * radius_interval, 3))));
    }
}

void termostatoAnderson(System &system, double freq, double dt, double sigma_velocities) {
    int n_particles_to_reset = floor(system.N_particles * freq * dt);

    for (size_t i = 0; i < n_particles_to_reset; i++) {
        int index = rand() / (RAND_MAX + 1.) * system.N_particles;
        system.particles[index].vel.x = getRNDVelocity(sigma_velocities);
        system.particles[index].vel.y = getRNDVelocity(sigma_velocities);
        system.particles[index].vel.z = getRNDVelocity(sigma_velocities);
    }
}

void termostatoAMoltiplicatore(System &system, double aim_temperature) {
    double current_temperature = TQ_Temperature(system);
    double lambda = sqrt(aim_temperature / current_temperature);
    for (size_t i = 0; i < system.N_particles; i++) {
        system.particles[i].vel.x *= lambda;
        system.particles[i].vel.y *= lambda;
        system.particles[i].vel.z *= lambda;
    }
}

void gasSimulation() {
    //------------------------------------
    const int n_init = 13;
    InitialCondition init_conditions[n_init] = {
        {.densita = 1.2, .temperatura = 2, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_solido.dat", .file_name_thermo = "./data/thermo_solido.dat"},
        {.densita = 0.01, .temperatura = 1, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_gas.dat", .file_name_thermo = "./data/thermo_gas.dat"},
        {.densita = 0.8, .temperatura = 1.8, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_liquido.dat", .file_name_thermo = "./data/thermo_liquido.dat"},
        {.densita = 0.7, .temperatura = 1.7, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.6, .temperatura = 1.3, .stampa = 0, .file_name_thermo = ""},
        {.densita = 1, .temperatura = 1.3, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.75, .temperatura = 1.8, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.2, .temperatura = 0.6, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.4, .temperatura = 0.65, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.001, .temperatura = 1, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.5, .temperatura = 0.9, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.3, .temperatura = 0.57, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.25, .temperatura = 0.55, .stampa = 0, .file_name_thermo = ""}};
    // -----------------------------------
    int N_PARTICLES = pow(N_CELLS, 3) * M;
    double t0 = 0.;    // Time zero simulation
    double tf = 6.;    // Time final simulation
    double dt = 1e-2;  // Time interval
    double N_time_steps = (tf - t0) / dt;
    //
    FILE *p_t_ro = fopen("./data/pressione_temperatura_rho.dat", "w");
    // SYSTEM INITIALIZATION
    System system;
    system.t = t0;
    system.N_particles = N_PARTICLES;

    for (size_t i = 0; i < n_init; i++) {
        srand(SEED);
        FILE *distribuzioneRadiale_file = fopen(init_conditions[i].file_name_g, "w");
        FILE *termo_data_file = fopen(init_conditions[i].file_name_thermo, "w");
        //  Struttura sistema
        system.densita = init_conditions[i].densita;
        double volume = (system.N_particles / system.densita);
        double L = cbrt(volume);
        system.volume = volume;
        system.L = L;
        system.particles = (Particle *)malloc(sizeof(Particle) * system.N_particles);

        // Temperature
        double sigma_velocities = sqrt(init_conditions[i].temperatura);

        for (size_t j = 0; j < system.N_particles; j++) {
            system.particles[j].mass = MASS;
        }

        printf("--------------------------------\n");
        printf("Velocities inizialized system %d/%d\n", i + 1, n_init);
        initVelocities(system, sigma_velocities);
        printf("Positions inizialized system %d/%d\n", i + 1, n_init);
        initPositions(system);
        system.forcesWork = (double *)calloc(system.N_particles, sizeof(double));
        system.forces = getForcesLennarJones(system, NULL);
        printf("Avvio simulazione con %d particelle\n", system.N_particles);
        // Initialize thermodinamic variables arrays
        double *temperature_array = (double *)malloc(sizeof(double) * N_time_steps);
        double *pressure_array = (double *)malloc(sizeof(double) * N_time_steps);
        double *energy_array = (double *)malloc(sizeof(double) * N_time_steps);

        // System evolution

        for (size_t j = 0; j < N_time_steps; j++) {
            temperature_array[j] = TQ_Temperature(system);
            pressure_array[j] = TQ_Pressure(system);
            energy_array[j] = system.kinetic_en + system.pot_en;
            verletPropagator(system, dt, getForcesLennarJones, NULL);

            // termostatoAnderson(system, FREQ, dt, sigma_velocities);
            //  termostatoAMoltiplicatore(system, init_conditions[i].temperatura);
            //     printf("Energy: %f\n", system.kinetic_en + system.pot_en);
        }

        int f_step = round(3 / dt);
        // Stampa termodinamica
        if (init_conditions[i].stampa) {
            printf("Calocolo distribuzione radiale.\n");
            distribuzioneRadiale(system, distribuzioneRadiale_file);
            for (size_t j = 0; j < N_time_steps; j++) {
                fprintf(termo_data_file, "%lf %lf %lf %lf\n", t0 + j * dt, temperature_array[j], pressure_array[j], energy_array[j]);
            }
        }

        printf("Temperatura: %f +- %f\n", mean_array(temperature_array, N_time_steps, f_step), sqrt(var_array(temperature_array, N_time_steps, f_step)));
        printf("Compressibilità: %f +- %f\n", mean_array(pressure_array, N_time_steps, f_step), sqrt(var_array(pressure_array, N_time_steps, f_step)));
        printf("Enengia: %f +- %f\n", mean_array(energy_array, N_time_steps, f_step), sqrt(var_array(energy_array, N_time_steps, f_step)));

        fprintf(p_t_ro, "%lf, %lf, %lf\n", system.densita, mean_array(temperature_array, N_time_steps, f_step), mean_array(pressure_array, N_time_steps, f_step));

        free(temperature_array);
        free(pressure_array);
        free(energy_array);
    }
}

void armonicOscillator() {
    double t0 = 0.;
    double tf = 1.;
    double dt = 1e-3;
    int N_steps = (tf - t0) / dt;

    FILE *armonicOscillator_verlet_file = fopen("./data/armonic_oscillator_velret.dat", "w");

    System system;
    system.N_particles = 1;
    system.particles = (Particle *)malloc(sizeof(Particle));
    system.particles[0].pos = Vec3{.x = 0, .y = 0, .z = 0};
    system.particles[0].vel = Vec3{.x = 1, .y = 0, .z = 0};
    system.particles[0].mass = 1;
    double k[] = {1, 1, 1};
    system.forces = getForcesOscillatore(system, k);
    for (size_t i = 0; i < N_steps; i++) {
        verletPropagator(system, dt, getForcesOscillatore, k);
        fprintf(armonicOscillator_verlet_file, "%lf %lf %lf %lf\n", system.t, system.particles[0].pos.x, system.particles[0].vel.x);
    }
}

int main(int argc, char const *argv[]) {
    gasSimulation();
    //  armonicOscillator();
}
