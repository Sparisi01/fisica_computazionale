#include <math.h>

#include <fstream>
#include <iostream>

struct vec3
{
    double t;
    double x;
    double y;
};

struct vecMR
{
    double R;
    double M;
};

using namespace std;

int euleroEsplicito(vec3 *q, double h, vec3 (*F)(vec3, double *), double *args)
{
    *q = vec3{
        q->t + h,
        q->x + h * F(*q, args).x,
        q->y + h * F(*q, args).y,
    };
    return 1;
}

int rk4(vec3 *q, double h, vec3 (*F)(vec3, double *), double *args)
{
    double k1, k2, k3, k4, l1, l2, l3, l4;
    vec3 tempF;

    k1 = h * F(*q, args).x;
    l1 = h * F(*q, args).y;

    tempF = F(vec3{q->t + h / 2, q->x + k1 / 2, q->y + l1 / 2}, args);
    k2 = h * tempF.x;
    l2 = h * tempF.y;

    tempF = F(vec3{q->t + h / 2, q->x + k2 / 2, q->y + l2 / 2}, args);
    k3 = h * tempF.x;
    l3 = h * tempF.y;

    tempF = F(vec3{q->t + h, q->x + k3, q->y + l3}, args);
    k4 = h * tempF.x;
    l4 = h * tempF.y;

    *q = vec3{
        q->t + h,
        q->x + (k1 + 2 * (k2 + k3) + k4) / 6.,
        q->y + (l1 + 2 * (l2 + l3) + l4) / 6.,
    };

    return 1;
}

double density(double P, double k, double gamma)
{
    return pow(P / (k * (gamma - 1)), 1 / gamma);
}

vec3 F_NR(vec3 q, double *parameters)
{
    return vec3{
        // r = t, x = m, y = p
        q.t,
        pow(q.t, 2) * density(q.y, parameters[0], parameters[1]),
        -q.x * density(q.y, parameters[0], parameters[1]) * pow(q.t, -2),
    };
}

vec3 F_R(vec3 q, double *parameters)
{
    double current_density = density(q.y, parameters[0], parameters[1]);
    double energy = current_density + parameters[0] * pow(current_density, parameters[1]);
    return vec3{
        q.t,
        pow(q.t, 2) * energy,
        -(q.y + energy) * (q.x + pow(q.t, 3) * q.y) / (pow(q.t, 2) - 2 * q.t * q.x),
    };
}

#define PI 3.141592653
// #define PI 3.14
#define C 299792458
#define R_max 100.
#define R0 (1 / sqrt(4 * PI * 197.327 * 6.67259 * 0.16 * 938.565 * 1e-45) * 1e-18) // R0 in km
// #define M0 (1 / sqrt(4 * PI * 0.16 * 938.565 * pow(197.327 * 6.67259 * 1e-45, 3)))  // M0 in Mev / c^2
#define M0 (R0 / (197.327 * 6.67259 * 1e-63 * 4 * PI))

vecMR find_M_R_mk4(double h, int N_steps, vec3 *q, vec3 (*F)(vec3, double *), double *par_star)
{
    double last_r = INFINITY;
    double last_m = INFINITY;
    for (size_t j = 0; j < N_steps; j++)
    {
        rk4(q, h, F, par_star);
        if (isnan(q->y) || q->y <= 0) break;
        last_r = q->t;
        last_m = q->x;
    }

    return vecMR{.R = last_r, .M = last_m};
}

vecMR find_M_R_eulero(double h, int N_steps, vec3 *q, vec3 (*F)(vec3, double *), double *par_star)
{
    double last_r = INFINITY;
    double last_m = INFINITY;
    for (size_t j = 0; j < N_steps; j++)
    {
        euleroEsplicito(q, h, F_NR, par_star);
        if (isnan(q->y) || q->y <= 0) break;
        last_r = q->t;
        last_m = q->x;
    }

    return vecMR{.R = last_r, .M = last_m};
}

int convergenza(double precisione = 1E-6)
{
    ofstream convergenza_rk4("./data/convergenza_rk4.dat");
    ofstream convergenza_eulero("./data/convergenza_eulero.dat");

    double P_center;

    double par_star_1[] = {/*k*/ 0.05, /*Gamma*/ 5. / 3.};
    double par_star_2[] = {0.1, 4. / 3.};
    double par_star_3[] = {0.01, 2.54};

    double N_h = 100;
    double N_k = 100;

    /*     double h = R_max / N_h;
        double k = R_max / N_k; */
    double h = 0.5;
    double k = 0.5;

    bool mk4_end = false;
    bool eulero_end = false;

    vecMR star_MR_mk4_last = {.R = INFINITY, .M = INFINITY};
    vecMR star_MR_eulero_last = {.R = INFINITY, .M = INFINITY};

    P_center = pow(2, -10);
    while (!eulero_end || !mk4_end)
    {
        vec3 initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        vec3 *q = &initialCondition;
        if (!mk4_end)
        {
            vecMR star_MR_mk4 = find_M_R_mk4(h /= 2, N_h *= 2, q, F_NR, par_star_1);
            double diff = sqrt(pow(star_MR_mk4.M - star_MR_mk4_last.M, 2) + pow(star_MR_mk4.R - star_MR_mk4_last.R, 2));
            convergenza_rk4 << h << " " << diff << endl;
            if (diff < precisione)
            {
                std::cout << "Step mk4: " << h << endl;
                mk4_end = true;
            }
            star_MR_mk4_last.M = star_MR_mk4.M;
            star_MR_mk4_last.R = star_MR_mk4.R;
        }
        if (!eulero_end)
        {
            initialCondition = {.t = 0.001, .x = 0, .y = P_center};
            vecMR star_MR_eulero = find_M_R_eulero(k /= 2, N_k *= 2, q, F_NR, par_star_1);
            // printf("%f", star_MR_eulero.M);
            //  double diff = sqrt(pow(star_MR_eulero.M - star_MR_eulero_last.M, 2) + pow(star_MR_eulero.R - star_MR_eulero_last.R, 2));
            double diff = fabs(star_MR_eulero.M - star_MR_eulero_last.M) + fabs(star_MR_eulero.R - star_MR_eulero_last.R);
            // double diffM = star_MR_eulero.M - star_MR_eulero_last.M
            convergenza_eulero << k << " " << diff << endl;
            if (diff < precisione)
            {
                std::cout << "Step eulero: " << k << endl;
                eulero_end = true;
            }
            star_MR_eulero_last.M = star_MR_eulero.M;
            star_MR_eulero_last.R = star_MR_eulero.R;
        }
    }
    return 1;
}

void graficoDavide()
{
    double h0 = 1e-6;

    FILE *file_errore_rk = fopen("./data/errore_rk.dat", "w");
    double par_star[] = {0.1, 4. / 3.};
    double h = h0;
    vec3 initialCondition = {.t = 0.001, .x = 0, .y = pow(2, -10)};
    vecMR v_MR_NR_0 = find_M_R_mk4(1e-9, R_max / 1e-9, &initialCondition, F_NR, par_star);
    for (size_t i = 0; i < 250; i++)
    {
        h *= 1.05;
        vec3 initialCondition = {.t = 0.001, .x = 0, .y = pow(2, -10)};
        vecMR v_MR_NR = find_M_R_mk4(h, R_max / h, &initialCondition, F_NR, par_star);
        double deltaRadius = fabs(v_MR_NR.R - v_MR_NR_0.R);
        fprintf(file_errore_rk, "%10.5E %10.5E %10.5E \n", h / h0, h, deltaRadius);
    }

    FILE *file_errore_eu = fopen("./data/errore_eu.dat", "w");
    h = h0;
    initialCondition = {.t = 0.001, .x = 0, .y = pow(2, -10)};
    v_MR_NR_0 = find_M_R_eulero(h, R_max / h, &initialCondition, F_NR, par_star);
    for (size_t i = 0; i < 250; i++)
    {
        h *= 1.05;
        vec3 initialCondition = {.t = 0.001, .x = 0, .y = pow(2, -10)};
        vecMR v_MR_NR = find_M_R_eulero(1e-9, R_max / 1e-9, &initialCondition, F_NR, par_star);
        double deltaRadius = fabs(v_MR_NR.R - v_MR_NR_0.R);
        fprintf(file_errore_eu, "%10.5E %10.5E %10.5E\n", h / h0, h, deltaRadius);
    }
}

int soluzioni(const char *filename, double P_center, vec3 (*F)(vec3, double *), double *parStar)
{
    ofstream file(filename);
    int N_steps = 10000;
    double h = R_max / N_steps;

    vec3 initialCondition = {.t = 0.001, .x = 0, .y = P_center};
    vec3 *q = &initialCondition;

    for (size_t j = 0; j < N_steps; j++)
    {
        rk4(q, h, F, parStar);
        if (isnan(q->y) || q->y <= 0) break;
        file << q->t << " "
             << " " << q->x << " " << q->y << endl;
    }
    return 1;
}

int main(int argc, char const *argv[])
{
    double M0_ = 13.6558;
    double R0_ = 20.0615;

    // Convergenza
    // convergenza(1E-6);
    graficoDavide();
    /* double par_star_1[] = {0.05, 5. / 3.};
    double par_star_2[] = {0.1, 4. / 3.};
    double par_star_3[] = {0.01, 2.54};

    // Draw soluzioni eq differenziale
    soluzioni("./data/sol_stella_1_NR.dat", 1, F_NR, par_star_1);
    soluzioni("./data/sol_stella_2_NR.dat", 1, F_NR, par_star_2);
    soluzioni("./data/sol_stella_3_NR.dat", 1, F_NR, par_star_3);

    soluzioni("./data/sol_stella_1_R.dat", 1, F_R, par_star_1);
    soluzioni("./data/sol_stella_2_R.dat", 1, F_R, par_star_2);
    soluzioni("./data/sol_stella_3_R.dat", 1, F_R, par_star_3);

    // M R al variare della pressione
    ofstream P_M_R_1_file("./data/P_M_R_1.dat");
    ofstream P_M_R_2_file("./data/P_M_R_2.dat");
    ofstream P_M_R_3_file("./data/P_M_R_3.dat");
    P_M_R_1_file << "#"
                 << "P centrale - massa NR - Raggio NR - Massa R - Raggio R" << endl;
    P_M_R_2_file << "#"
                 << "P centrale - massa NR - Raggio NR - Massa R - Raggio R" << endl;
    P_M_R_3_file << "#"
                 << "P centrale - massa NR - Raggio NR - Massa R - Raggio R" << endl;
    vec3 initialCondition;
    vecMR v_MR_R;
    vecMR v_MR_NR;
    for (double P_center = pow(2, -30); P_center < pow(2, 30); P_center *= 2)
    {
        // Stella 1
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_NR = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_NR, par_star_1);
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_R = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_R, par_star_1);

        P_M_R_1_file
            << P_center << " " << v_MR_NR.M * M0_ << " " << v_MR_NR.R * R0_ << " " << v_MR_R.M * M0_ << " " << v_MR_R.R * R0_ << endl;
        // Stella 2
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_NR = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_NR, par_star_2);
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_R = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_R, par_star_2);

        P_M_R_2_file
            << P_center << " " << v_MR_NR.M * M0_ << " " << v_MR_NR.R * R0_ << " " << v_MR_R.M * M0_ << " " << v_MR_R.R * R0_ << endl;
        // Stella 3
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_NR = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_NR, par_star_3);
        initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        v_MR_R = find_M_R_mk4(1e-3, 1e3 * R_max, &initialCondition, F_R, par_star_3);

        P_M_R_3_file
            << P_center << " " << v_MR_NR.M * M0_ << " " << v_MR_NR.R * R0_ << " " << v_MR_R.M * M0_ << " " << v_MR_R.R * R0_ << endl;
    }

    std::cout << "Valore di R0: " << R0 << " [km]" << endl
              << "Valore di M0: " << M0 << " [MeV*c^-2]"
              << "oppure " << (M0 * pow(C, -2)) * (1.602 * 1E-13) << " [Kg]" << endl; */

    return 0;
}
