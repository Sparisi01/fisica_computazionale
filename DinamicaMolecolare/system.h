typedef struct System
{
    double t;
    int n_particles;
    double *x, *y, *z;
    double *vx, *vy, *vz;
    double *masses;
    double *forces;
    double *last_time_forces;
};
