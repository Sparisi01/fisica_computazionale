#include <math.h>

#include <fstream>
#include <iostream>

struct vec3 {
    double t;
    double x;
    double y;
};

struct vecMR {
    double R;
    double M;
};

using namespace std;

int euleroEsplicito(vec3* q, double h, vec3 (*F)(vec3, double*), double* args) {
    *q = vec3{
        q->t + h,
        q->x + h * F(*q, args).x,
        q->y + h * F(*q, args).y,
    };
    return 1;
}

int rk4(vec3* q, double h, vec3 (*F)(vec3, double*), double* args) {
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

double density(double P, double k, double gamma) {
    return pow(P / (k * (gamma - 1)), 1 / gamma);
}

vec3 F(vec3 q, double* parameters) {
    return vec3{
        q.t,
        pow(q.t, 2) * density(q.y, parameters[0], parameters[1]),
        -q.x * density(q.y, parameters[0], parameters[1]) * pow(q.t, -2),
    };
}

#define PI 3.141592653
#define C 299792458
#define R_max 10.
#define R0 (1 / sqrt(4 * PI * 197.327 * 6.67259 * 0.16 * 938.565 * 1e-45) * 1e-18)  // R0 in km
#define M0 (1 / sqrt(4 * PI * 0.16 * 938.565 * pow(197.327 * 6.67259 * 1e-45, 3)))  // M0 in Mev / c^2

vecMR find_M_R_mk4(double h, int N_steps, double P_center, vec3* q, double* par_star) {
    double last_r = INFINITY;
    double last_m = INFINITY;
    for (size_t j = 0; j < N_steps; j++) {
        last_r = q->t;
        last_m = q->x;
        rk4(q, h, F, par_star);
        if (isnan(q->y) || q->y <= 0) break;
    }

    return vecMR{.R = last_r, .M = last_m};
}

vecMR find_M_R_eulero(double h, int N_steps, double P_center, vec3* q, double* par_star) {
    double last_r = INFINITY;
    double last_m = INFINITY;
    for (size_t j = 0; j < N_steps; j++) {
        last_r = q->t;
        last_m = q->x;
        euleroEsplicito(q, h, F, par_star);
        if (isnan(q->y) || q->y <= 0) break;
    }

    return vecMR{.R = last_r, .M = last_m};
}

int convergenza(double precisione = 1E-6) {
    ofstream convergenza_rk4("convergenza_rk4.dat");
    ofstream convergenza_eulero("convergenza_eulero.dat");

    double P_center;

    double par_star_1[] = {/*k*/ 0.05, /*Gamma*/ 5. / 3.};
    double par_star_2[] = {0.1, 4. / 3.};
    double par_star_3[] = {0.01, 2.54};

    double N_h = 100;
    double N_k = 100;

    double h = R_max / N_h;
    double k = R_max / N_k;

    bool mk4_end = false;
    bool eulero_end = false;

    vecMR star_MR_mk4_last = {.R = INFINITY, .M = INFINITY};
    vecMR star_MR_eulero_last = {.R = INFINITY, .M = INFINITY};

    P_center = pow(2, -10);
    while (!eulero_end || !mk4_end) {
        vec3 initialCondition = {.t = 0.001, .x = 0, .y = P_center};
        vec3* q = &initialCondition;
        if (!mk4_end) {
            vecMR star_MR_mk4 = find_M_R_mk4(h /= 2, N_h *= 2, P_center, q, par_star_1);
            double diff = sqrt(pow(star_MR_mk4.M - star_MR_mk4_last.M, 2) + pow(star_MR_mk4.R - star_MR_mk4_last.R, 2));
            convergenza_rk4 << h << " " << diff << endl;
            if (diff < precisione) {
                std::cout << "Step mk4: " << h << endl;
                mk4_end = true;
            }
            star_MR_mk4_last.M = star_MR_mk4.M;
            star_MR_mk4_last.R = star_MR_mk4.R;
        }
        if (!eulero_end) {
            initialCondition = {.t = 0.001, .x = 0, .y = P_center};
            vecMR star_MR_eulero = find_M_R_eulero(k /= 2, N_k *= 2, P_center, q, par_star_1);
            double diff = sqrt(pow(star_MR_eulero.M - star_MR_eulero_last.M, 2) + pow(star_MR_eulero.R - star_MR_eulero_last.R, 2));
            convergenza_eulero << k << " " << diff << endl;
            if (diff < precisione) {
                std::cout << "Step eulero: " << k << endl;
                eulero_end = true;
            }
            star_MR_eulero_last.M = star_MR_eulero.M;
            star_MR_eulero_last.R = star_MR_eulero.R;
        }
    }
    return 1;
}

int soluzioni(const char* filename, double* parStar, double P_center) {
    ofstream file(filename);
    int N_steps = 10000;
    double h = R_max / N_steps;

    vec3 initialCondition = {.t = 0.001, .x = 0, .y = P_center};
    vec3* q = &initialCondition;

    for (size_t j = 0; j < N_steps; j++) {
        rk4(q, h, F, parStar);
        if (isnan(q->y) || q->y <= 0) break;
        file << q->t << " "
             << " " << q->x << " " << q->y << endl;
    }
    return 1;
}

int main(int argc, char const* argv[]) {
    // Convergenza
    //convergenza(1E-6);

    // Draw soluzioni
    double par_star_1[] = {/*k*/ 0.05, /*Gamma*/ 5. / 3.};
    double par_star_2[] = {0.1, 4. / 3.};
    double par_star_3[] = {0.01, 2.54};
    soluzioni("stella1.dat", par_star_1, 4);
    soluzioni("stella2.dat", par_star_2, 4);
    soluzioni("stella3.dat", par_star_3, 4);

    std::cout << "Valore di R0: " << R0 << " [km]" << endl
              << "Valore di M0: " << M0 << " [MeV*c^-2]"
              << "oppure " << (M0 * pow(C, -2)) * (1.602 * 1E-13) << " [Kg]" << endl;

    return 0;
}
