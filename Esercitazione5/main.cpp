#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <iostream>

using namespace std;

struct Vec3
{
    double x, y, z;
};
struct Particle
{
    Vec3 pos, vel; // Particle pos and vel vector
    double mass;   // Particle mass
};

struct System
{
    double t;            // Time
    int N_particles;     // Number of particles
    Particle *particles; // Vector of particle
    Vec3 *forces;        // Force[i] vector acting on particle[i]
    double *forcesWork;  // Sum of all internal force work ΣFᵢ(t)∙rᵢ(t)
    double kinetic_en;   // Kinetic_energy(t)
    double pot_en;       // U(t)
    double density;      // Particles density N_particles/(L*L*L)
    double L;            // Box lenght
};

struct InitialCondition
{
    double densita;
    double temperatura;
    int stampa;
    const char *file_name_g;
    const char *file_name_thermo;
};

struct Thermodinamics
{
    double mean_T, sigma_T;
    double mean_P, sigma_P;
    double mean_E, sigma_E;
};

// TEMP - RHO - M - N_CELLS - TIME_CUT - TEq
// 2 - 1.2 - 4 - 4 - 0.4 - 1.062724

#define PI 3.14159265359
#define SEED 3

// UNITÀ DI MISURA
#define SIGMA 1.
#define EPSILON 1.
#define MASS 1.
// STRUTTURA RETICOLO
#define M 4       // M=1 CC, M=2 BCC, M=4 FCC
#define N_CELLS 5 // Numero celle per lato
#define FREQ 0    // Frequenza termostato, 0 disattivato
// OUTPUT
#define PRINT_THERMO 1

void verletPropagator(System &system, double dt, Vec3 *(*F)(System &, double *), double *args)
{

    system.kinetic_en = 0;
    system.pot_en = 0;

    Vec3 *oldForces = system.forces;
    // Aggiorna posizioni di tutte le particelle
    for (size_t i = 0; i < system.N_particles; i++)
    {
        system.particles[i].pos.x += system.particles[i].vel.x * dt + 0.5 / system.particles[i].mass * oldForces[i].x * dt * dt;
        system.particles[i].pos.y += system.particles[i].vel.y * dt + 0.5 / system.particles[i].mass * oldForces[i].y * dt * dt;
        system.particles[i].pos.z += system.particles[i].vel.z * dt + 0.5 / system.particles[i].mass * oldForces[i].z * dt * dt;
    }
    // Incrementa il tempo del sistema
    system.t += dt;

    // Calcola forze agenti sulle particelle
    Vec3 *newForces = F(system, args);

    //  Aggiorna velocità
    for (size_t i = 0; i < system.N_particles; i++)
    {
        system.kinetic_en += system.particles[i].mass * (system.particles[i].vel.x * system.particles[i].vel.x) / 2.;
        system.kinetic_en += system.particles[i].mass * (system.particles[i].vel.y * system.particles[i].vel.y) / 2.;
        system.kinetic_en += system.particles[i].mass * (system.particles[i].vel.z * system.particles[i].vel.z) / 2.;

        system.particles[i].vel.x += 0.5 / system.particles[i].mass * (oldForces[i].x + newForces[i].x) * dt;
        system.particles[i].vel.y += 0.5 / system.particles[i].mass * (oldForces[i].y + newForces[i].y) * dt;
        system.particles[i].vel.z += 0.5 / system.particles[i].mass * (oldForces[i].z + newForces[i].z) * dt;
    }

    system.forces = newForces;
    free(oldForces);
}

double lennarJonesPotential(const System &system, double r)
{
    if (r >= system.L / 2)
    {
        return 0;
    }
    else
    {
        double VSHIFT = (4 * EPSILON * (pow(SIGMA / (system.L / 2), 12) - pow(SIGMA / (system.L / 2), 6)));
        return 4. * EPSILON * (pow(SIGMA / r, 12) - pow(SIGMA / r, 6)) - VSHIFT;
    }
}

Vec3 *getForcesLennarJones(System &system, double *args)
{
    struct Vec3 *forces = (Vec3 *)calloc(system.N_particles, sizeof(struct Vec3));

    for (size_t i = 0; i < system.N_particles; i++)
    {
        forces[i].x = 0;
        forces[i].y = 0;
        forces[i].z = 0;
    }

    // Reset lavori forze
    for (size_t i = 0; i < system.N_particles; i++)
        system.forcesWork[i] = 0;

    // Cicli for per ogni particella
    for (size_t i = 0; i < system.N_particles - 1; i++)
    {
        for (size_t j = i + 1; j < system.N_particles; j++)
        {
            Vec3 &p1 = system.particles[i].pos;
            Vec3 &p2 = system.particles[j].pos;
            double x = (p1.x - p2.x) - system.L * rint((p1.x - p2.x) / system.L);
            double y = (p1.y - p2.y) - system.L * rint((p1.y - p2.y) / system.L);
            double z = (p1.z - p2.z) - system.L * rint((p1.z - p2.z) / system.L);
            double r_ij_mod = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
            Vec3 r_dir_ij = Vec3{.x = x / r_ij_mod, .y = y / r_ij_mod, .z = z / r_ij_mod};
            // double F_magnitude = 24 * EPSILON * SIGMA / (r * r) * (2 * pow(SIGMA / r, 11) - pow(SIGMA / r, 5));

            // Use formula F = -dV/dr
            double h = 1e-6;
            double F_magnitude = -(lennarJonesPotential(system, r_ij_mod + h) - lennarJonesPotential(system, r_ij_mod - h)) / (2 * h);

            forces[i].x += r_dir_ij.x * F_magnitude;
            forces[i].y += r_dir_ij.y * F_magnitude;
            forces[i].z += r_dir_ij.z * F_magnitude;

            forces[j].x -= r_dir_ij.x * F_magnitude;
            forces[j].y -= r_dir_ij.y * F_magnitude;
            forces[j].z -= r_dir_ij.z * F_magnitude;

            // system.forcesWork[i] += force_ij * (r_dir_ij * r);
            // system.forcesWork[j] += force_ij * (r_dir_ij * r);

            system.forcesWork[i] += F_magnitude * r_ij_mod;
            system.pot_en += lennarJonesPotential(system, r_ij_mod);
        }
    }

    return forces;
}

double getRNDVelocity(double SIGMA_VELOCITIES)
{
    double x = rand() / (RAND_MAX + 1.);
    double y = rand() / (RAND_MAX + 1.);

    return SIGMA_VELOCITIES * sqrt(-2 * log(1 - x)) * cos(2 * PI * y);
}

// Initialize velocities
void initVelocities(System &system, double SIGMA_VELOCITIES)
{
    srand(SEED);
    for (size_t i = 0; i < system.N_particles; i += 1)
    {
        system.particles[i].vel.x = getRNDVelocity(SIGMA_VELOCITIES);
        system.particles[i].vel.y = getRNDVelocity(SIGMA_VELOCITIES);
        system.particles[i].vel.z = getRNDVelocity(SIGMA_VELOCITIES);
    }
}

// Initialize positions in a cubic lattice
void initPositions(System &system)
{
    Vec3 lattice_positions[M];
    switch (M)
    {
    case 1: // CC
        lattice_positions[0] = {.x = 0, .y = 0, .z = 0};
        break;
    case 2: // BCC
        lattice_positions[0] = {.x = 0, .y = 0, .z = 0};
        lattice_positions[1] = {.x = 0.5, .y = 0.5, .z = 0.5};
        break;
    case 4: // FCC
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
    for (size_t i = 0; i < N_CELLS; i++)
    {
        for (size_t j = 0; j < N_CELLS; j++)
        {
            for (size_t k = 0; k < N_CELLS; k++)
            {
                for (size_t w = 0; w < M; w++)
                {
                    system.particles[q].pos = Vec3{.x = (i + lattice_positions[w].x) * jump, .y = (j + lattice_positions[w].y) * jump, .z = (k + lattice_positions[w].z) * jump};
                    q++;
                }
            }
        }
    }
}

double mean_array(const double *array, int len, int first)
{
    double sum = 0;
    for (size_t i = first; i < len; i++)
    {
        sum += array[i];
    }

    return sum / (len - first);
}

double var_array(const double *array, int len, int first)
{
    double mean = mean_array(array, len, first);
    double sum = 0;
    for (size_t i = first; i < len; i++)
    {
        sum += pow(array[i] - mean, 2);
    }

    return sum / (len - first);
}

double TQ_Temperature(const System &system)
{
    double sum = 0;
    for (size_t i = 0; i < system.N_particles; i++)
    {
        sum += system.particles[i].mass * pow(system.particles[i].vel.x, 2);
        sum += system.particles[i].mass * pow(system.particles[i].vel.y, 2);
        sum += system.particles[i].mass * pow(system.particles[i].vel.z, 2);
    }
    return sum / (3 * system.N_particles);
}

double TQ_Pressure(const System &system)
{
    double W = 0;
    for (size_t i = 0; i < system.N_particles; i++)
    {
        W += system.forcesWork[i];
    }
    W /= system.N_particles;
    double T = TQ_Temperature(system);

    // return (T * system.density + W / 3); // Pressione
    return (1 + W / (3 * T)); // Compressibiità
}

void termostatoAnderson(System &system, double freq, double dt, double sigma_velocities)
{
    int n_particles_to_reset = floor(system.N_particles * freq * dt);

    for (size_t i = 0; i < n_particles_to_reset; i++)
    {
        int index = rand() / (RAND_MAX + 1.) * system.N_particles;
        system.particles[index].vel.x = getRNDVelocity(sigma_velocities);
        system.particles[index].vel.y = getRNDVelocity(sigma_velocities);
        system.particles[index].vel.z = getRNDVelocity(sigma_velocities);
    }
}

void termostatoAMoltiplicatore(System &system, double aim_temperature)
{
    double current_temperature = TQ_Temperature(system);
    double lambda = sqrt(aim_temperature / current_temperature);
    for (size_t i = 0; i < system.N_particles; i++)
    {
        system.particles[i].vel.x *= lambda;
        system.particles[i].vel.y *= lambda;
        system.particles[i].vel.z *= lambda;
    }
}

void writePositions(const System &system, FILE *file)
{
    for (size_t i = 0; i < system.N_particles; i++)
    {
        fprintf(file, "%lf %lf %lf\n", system.particles[i].pos.x, system.particles[i].pos.y, system.particles[i].pos.z);
    }
}

Thermodinamics gasSimulation(InitialCondition init_condition)
{

    // -----------------------------------
    int N_PARTICLES = pow(N_CELLS, 3) * M;
    const double t0 = 0.;   // Time zero simulation
    const double tf = 6.;   // Time final simulation
    const double dt = 1e-2; // Time interval
    const double N_time_steps = (tf - t0) / dt;

    // SYSTEM INITIALIZATION
    System system;
    system.N_particles = N_PARTICLES;

    FILE *distribuzioneRadiale_file = fopen(init_condition.file_name_g, "w");
    FILE *termo_data_file = fopen(init_condition.file_name_thermo, "w");

    system.t = t0;
    system.density = init_condition.densita;
    double volume = (system.N_particles / system.density); // Cell volume
    system.L = cbrt(volume);
    system.particles = (Particle *)malloc(sizeof(Particle) * system.N_particles);
    if (!system.particles)
    {
        perror("Error allocating particles work array");
        exit(EXIT_FAILURE);
    }

    double sigma_velocities = sqrt(init_condition.temperatura);

    // Init masses
    for (size_t j = 0; j < system.N_particles; j++)
    {
        system.particles[j].mass = MASS;
    }

    initVelocities(system, sigma_velocities);
    initPositions(system);

    // FILE *file_init_position = fopen("./data/posizione_init.dat", "w");
    // writePositions(system, file_init_position);

    system.forcesWork = (double *)calloc(system.N_particles, sizeof(double));
    if (!system.forcesWork)
    {
        perror("Error allocating force work array");
        exit(EXIT_FAILURE);
    }

    system.forces = getForcesLennarJones(system, NULL);

    printf("N Particles: %d\n", system.N_particles);
    printf("Box length: %lf\n", system.L);

    // Initialize thermodinamic variables arrays
    double *temperature_array = (double *)malloc(sizeof(double) * N_time_steps);
    double *pressure_array = (double *)malloc(sizeof(double) * N_time_steps);
    double *energy_array = (double *)malloc(sizeof(double) * N_time_steps);
    if (!temperature_array || !pressure_array || !energy_array)
    {
        perror("Error allocating thermo array");
        exit(EXIT_FAILURE);
    }

    int f_step = round(3 / dt);   // Primo step da cui iniziare a fare statistica
    double max_radius = system.L; // Raggio massimo g(r)
    double N_intervals = 600;     // N bins
    double radius_interval = max_radius / N_intervals;

    int *counting_array = (int *)calloc(sizeof(int), N_intervals); // Array di bin per tenere traccia della densità radiale
    if (!counting_array)
    {
        perror("Error allocating counting array");
        exit(EXIT_FAILURE);
    }
    int n_g_rad_done = 0; // Numero di densità radiali calcolare, usato per la normalizzazione della g

    for (size_t j = 0; j < N_time_steps; j++)
    {
        // Thermodinamic at time t
        temperature_array[j] = TQ_Temperature(system);
        pressure_array[j] = TQ_Pressure(system);
        energy_array[j] = system.kinetic_en + system.pot_en;

        // Bring particles back in the (0,0,0) box
        for (size_t i = 0; i < system.N_particles; i++)
        {
            system.particles[i].pos.x = system.particles[i].pos.x - system.L * rint(system.particles[i].pos.x / system.L);
            system.particles[i].pos.y = system.particles[i].pos.y - system.L * rint(system.particles[i].pos.y / system.L);
            system.particles[i].pos.z = system.particles[i].pos.z - system.L * rint(system.particles[i].pos.z / system.L);
        }

        verletPropagator(system, dt, getForcesLennarJones, NULL);

        // termostatoAnderson(system, FREQ, dt, sigma_velocities);
        // termostatoAMoltiplicatore(system, init_conditions[i].temperatura);
        // printf("Energy: %f\n", system.kinetic_en + system.pot_en);

        // Compute g(r,t), only after thermalization
        if (j > f_step && init_condition.stampa)
        {
            n_g_rad_done++;
            for (size_t i = 0; i < system.N_particles; i++)
            {
                Vec3 &p1 = system.particles[i].pos;
                for (size_t k = 0; k < system.N_particles; k++)
                {
                    if (k == i) continue;
                    double x = (p1.x - system.particles[k].pos.x) - system.L * rint((p1.x - system.particles[k].pos.x) / system.L);
                    double y = (p1.y - system.particles[k].pos.y) - system.L * rint((p1.y - system.particles[k].pos.y) / system.L);
                    double z = (p1.z - system.particles[k].pos.z) - system.L * rint((p1.z - system.particles[k].pos.z) / system.L);
                    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
                    int n_radius_jump = (int)round(r / radius_interval);
                    counting_array[n_radius_jump]++;
                }
            }
        }
    }

    // Stampa termodinamica
    if (init_condition.stampa)
    {
        for (size_t j = 0; j < N_intervals; j++)
        {

            fprintf(distribuzioneRadiale_file, "%f %f %f\n", j * radius_interval, j * radius_interval * 2 / system.L, (double)(counting_array[j]) / (system.density * system.N_particles * n_g_rad_done * (4 * PI / 3) * (pow((j + 1) * radius_interval, 3) - pow(j * radius_interval, 3))));
        }
        for (size_t j = 0; j < N_time_steps; j++)
        {
            fprintf(termo_data_file, "%lf %lf %lf %lf\n", t0 + j * dt, temperature_array[j], pressure_array[j], energy_array[j]);
        }
    }

    Thermodinamics thermo = {
        .mean_T = mean_array(temperature_array, N_time_steps, f_step),
        .sigma_T = sqrt(var_array(temperature_array, N_time_steps, f_step)),
        .mean_P = mean_array(pressure_array, N_time_steps, f_step),
        .sigma_P = sqrt(var_array(pressure_array, N_time_steps, f_step)),
        .mean_E = mean_array(energy_array, N_time_steps, f_step),
        .sigma_E = sqrt(var_array(energy_array, N_time_steps, f_step)),
    };

    // FILE *file_end_position = fopen("./data/posizione_end.dat", "w");
    // writePositions(system, file_end_position);
    // fclose(file_end_position);

    free(temperature_array);
    free(pressure_array);
    free(energy_array);
    free(system.particles);
    free(system.forcesWork);
    free(counting_array);
    fclose(termo_data_file);
    fclose(distribuzioneRadiale_file);

    return thermo;
}

int main(int argc, char const *argv[])
{

    //------------------------------------
    const int n_init = 3;
    InitialCondition init_conditions[n_init] = {
        {.densita = 1.2, .temperatura = 2.1, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_solido.dat", .file_name_thermo = "./data/thermo_solido.dat"},
        {.densita = 0.01, .temperatura = 1.1, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_gas.dat", .file_name_thermo = "./data/thermo_gas.dat"},
        {.densita = 0.8, .temperatura = 1.9, .stampa = 1, .file_name_g = "./data/distribuzione_radiale_liquido.dat", .file_name_thermo = "./data/thermo_liquido.dat"},
        /* {.densita = 0.7, .temperatura = 1.5, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.6, .temperatura = 1.15, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.1, .temperatura = 0.7, .stampa = 0, .file_name_thermo = ""},
        {.densita = 1, .temperatura = 1.3, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.75, .temperatura = 1.8, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.2, .temperatura = 0.45, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.4, .temperatura = 0.65, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.001, .temperatura = 1, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.5, .temperatura = 0.9, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.3, .temperatura = 0.57, .stampa = 0, .file_name_thermo = ""},
        {.densita = 0.25, .temperatura = 0.55, .stampa = 0, .file_name_thermo = ""} */
    };

    FILE *p_t_ro = fopen("./data/pressione_temperatura_rho.dat", "w");

    // LOOP over all init conditions and save themodinamics observables
    for (size_t i = 0; i < n_init; i++)
    {
        printf("-----------------------\n");
        printf("Starting simulation: %d/%d\n", i + 1, n_init);

        Thermodinamics thermo = gasSimulation(init_conditions[i]);

        printf("-----------------------\n");
        printf("Temperature:%f +- %f\n", thermo.mean_T, thermo.sigma_T);
        printf("Compressibility: %f +- %f\n", thermo.mean_P, thermo.sigma_P);
        printf("Energy: %f +- %f\n", thermo.mean_E, thermo.sigma_E);

        fprintf(p_t_ro, "%lf %lf %lf\n", init_conditions[i].densita, thermo.mean_T, thermo.mean_P);
    }

    // armonicOscillator();
}

/* Vec3 *getForcesOscillatore(System &system, double *args)
{
    struct Vec3 *forces = (Vec3 *)malloc(sizeof(struct Vec3) * system.N_particles);
    double k1 = args[0];
    double k2 = args[1];
    double k3 = args[2];
    for (size_t i = 0; i < system.N_particles; i++)
    {
        forces[i].x = -k1 * system.particles[i].pos.x;
        forces[i].y = -k2 * system.particles[i].pos.y;
        forces[i].z = -k3 * system.particles[i].pos.z;
    }
    return forces;
} */

/* void armonicOscillator()
{
    double t0 = 0.;
    double tf = 10.;
    double dt = 1e-3;
    int N_steps = (tf - t0) / dt;

    FILE *armonicOscillator_verlet_file = fopen("./data/armonic_oscillator_velret.dat", "w");

    System system;
    system.N_particles = 1;
    system.particles = (Particle *)malloc(sizeof(Particle));
    system.particles[0].pos = Vec3{.x = 1, .y = 0, .z = 0};
    system.particles[0].vel = Vec3{.x = 0, .y = 0, .z = 0};
    system.particles[0].mass = 1;
    double k[] = {1, 1, 1};
    system.forces = getForcesOscillatore(system, k);
    for (size_t i = 0; i < N_steps; i++)
    {
        verletPropagator(system, dt, getForcesOscillatore, k);
        fprintf(armonicOscillator_verlet_file, "%lf %lf %lf %lf\n", system.t, system.particles[0].pos.x, system.particles[0].vel.x);
    }
} */

/* void distribuzioneRadiale(const System &system, FILE *file)
{
    double max_radius = system.L / 2;
    double N_intervals = 300;
    double radius_interval = max_radius / N_intervals;
    int *counting_array = (int *)calloc(sizeof(int), N_intervals);
    // -----------------------------------
    for (size_t j = 0; j < system.N_particles; j++)
    {
        Vec3 &p1 = system.particles[j].pos;
        for (size_t i = 0; i < system.N_particles; i++)
        {
            if (i == j) continue;

            double x = (p1.x - system.particles[i].pos.x) - system.L * rint((p1.x - system.particles[i].pos.x) / system.L);
            double y = (p1.y - system.particles[i].pos.y) - system.L * rint((p1.y - system.particles[i].pos.y) / system.L);
            double z = (p1.z - system.particles[i].pos.z) - system.L * rint((p1.z - system.particles[i].pos.z) / system.L);
            double r_ij = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
            // printf("%lf", r / radius_interval);
            int n_radius_jump = round(r_ij / radius_interval);
            counting_array[n_radius_jump]++;
        }
    }
    for (size_t i = 0; i < N_intervals; i++)
    {
        // fprintf(file, "%f %d\n", i * radius_interval * 2 / system.L, counting_array[i]);
        fprintf(file, "%f %f\n", i * radius_interval * 2 / system.L, (double)(counting_array[i]) / (system.densita * system.N_particles * (4 * PI / 3) * (pow((i + 1) * radius_interval, 3) - pow(i * radius_interval, 3))));
    }
} */