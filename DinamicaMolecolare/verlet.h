#include "system.h";

void verletPropagator(System *system, double dt, double *(*F)(System, double *), double *args)
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
    // printf("%d ms \n", duration_cast<milliseconds>(tend - tin));
    //  Aggiorna velocità
    for (size_t i = 0; i < system.N_particles; i++)
    {
        system.kinetic_en += (system.particles[i].vel * system.particles[i].vel) / (2 * system.particles[i].mass);

        system.particles[i].vel.x += 0.5 / system.particles[i].mass * (oldForces[i].x + newForces[i].x) * dt;
        system.particles[i].vel.y += 0.5 / system.particles[i].mass * (oldForces[i].y + newForces[i].y) * dt;
        system.particles[i].vel.z += 0.5 / system.particles[i].mass * (oldForces[i].z + newForces[i].z) * dt;
    }

    system.forces = newForces;
    free(oldForces);
    // oldForces = NULL;
}