#include <math.h>

#include <fstream>
#include <iostream>

struct vec3 {
    double t;
    double x;
    double y;
};

using namespace std;

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
#define N 10000
#define R_max 10.
#define R0 (1 / sqrt(4 * PI * 197.327 * 6.67259 * 0.16 * 938.565 * 1e-45) * 1e-18)  // R0 in km
#define M0 (1 / sqrt(4 * PI * 0.16 * 938.565 * pow(197.327 * 6.67259 * 1e-45, 3)))  // M0 in Mev / c^2

int main(int argc, char const* argv[]) {
    ofstream radius_mass_data1("dati_stelle1.dat");
    ofstream radius_mass_data2("dati_stelle2.dat");
    ofstream radius_mass_data3("dati_stelle3.dat");

    double h = R_max / N;
    double P_center = pow(2, -10);

    double par_star_1[] = {/*k*/ 0.05, /*Gamma*/ 5. / 3.};
    double par_star_2[] = {0.1, 4. / 3.};
    double par_star_3[] = {0.01, 2.54};

    int N_stars = 20;

    double last_r = INFINITY;
    double last_m = INFINITY;

    // Stelle tipo 1
    for (size_t i = 0; i < N_stars; i++) {
        vec3 initialCondition = {/*Raggio zero*/ 0.01, /*Massa a raggio 0*/ 0, /*Pressione a raggio 0*/ P_center *= 2};
        vec3* q = &initialCondition;
        for (size_t j = 0; j < N; j++) {
            last_r = q->t;
            last_m = q->x;
            rk4(q, h, F, par_star_1);
            if (isnan(q->y) || q->y <= 0) break;
        }
        radius_mass_data1 << last_r << " " << last_m << endl;
    }

    // Stelle tipo 2
    P_center = pow(2, -10);
    for (size_t i = 0; i < N_stars; i++) {
        vec3 initialCondition = {/*Raggio zero*/ 0.01, /*Massa a raggio 0*/ 0, /*Pressione a raggio 0*/ P_center *= 2};
        vec3* q = &initialCondition;
        for (size_t j = 0; j < N; j++) {
            last_r = q->t;
            last_m = q->x;
            rk4(q, h, F, par_star_2);
            if (isnan(q->y) || q->y <= 0) break;
        }
        radius_mass_data2 << last_r << " " << last_m << endl;
    }

    // Stelle tipo 3
    P_center = pow(2, -17);
    for (size_t i = 0; i < N_stars; i++) {
        vec3 initialCondition = {/*Raggio zero*/ 0.01, /*Massa a raggio 0*/ 0, /*Pressione a raggio 0*/ P_center *= 2};
        vec3* q = &initialCondition;
        for (size_t j = 0; j < N; j++) {
            last_r = q->t;
            last_m = q->x;
            rk4(q, h, F, par_star_3);
            if (isnan(q->y) || q->y <= 0) break;
        }
        radius_mass_data3 << last_r << " " << last_m << endl;
    }

    cout << "Valore di R0: " << R0 << " [km]" << endl
         << "Valore di M0: " << M0 << " [MeV*c^-2]" << endl
         << "Valore di M0: " << (M0 / 299792458 / 299792458) * (1.602 * 1e-13) << " [Kg]" << endl;

    return 0;
}
