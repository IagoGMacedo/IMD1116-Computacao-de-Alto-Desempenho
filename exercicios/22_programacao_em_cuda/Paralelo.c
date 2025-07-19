#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <float.h>

#define Lx 1.0
#define Ly 1.0
#define Lz 1.0
#define Nx 241  // Igual em ambos!
#define Ny 241  // Igual em ambos!
#define Nz 241  // Igual em ambos!
#define NU 0.1
#define DT 0.00002
#define T 0.01
#define SIGMA 0.1

double elapsed(struct timeval t0, struct timeval t1) {
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1e6;
}

typedef struct {
    double mass;
    double max;
    double min;
    double l2_squared;
} Stats;

void compute_stats(double* u, Stats* stats) {
    double local_mass = 0.0;
    double local_max = -DBL_MAX;
    double local_min = DBL_MAX;
    double local_l2 = 0.0;

    #pragma omp parallel for reduction(+:local_mass,local_l2) \
                             reduction(max:local_max) \
                             reduction(min:local_min) \
                             collapse(3) schedule(static)
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                double val = u[i*Ny*Nz + j*Nz + k];
                local_mass += val;
                local_max = fmax(local_max, val);
                local_min = fmin(local_min, val);
                local_l2 += val * val;
            }
        }
    }

    stats->mass = local_mass;
    stats->max = local_max;
    stats->min = local_min;
    stats->l2_squared = local_l2;
}

void print_stats(Stats* stats, int step) {
    printf("Passo %d:\n", step);
    printf("  Massa total: %.6f\n", stats->mass);
    printf("  Valor máximo: %.6f\n", stats->max);
    printf("  Valor mínimo: %.6f\n", stats->min);
    printf("  Norma L2: %.6f\n\n", sqrt(stats->l2_squared));
}

int main() {
    int nt = (int)(T / DT);
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);
    double dz = Lz / (Nz - 1);

    double *u = malloc(Nx * Ny * Nz * sizeof(double));
    double *u_new = malloc(Nx * Ny * Nz * sizeof(double));
    double *x = malloc(Nx * sizeof(double));
    double *y = malloc(Ny * sizeof(double));
    double *z = malloc(Nz * sizeof(double));
    Stats stats;

    // Inicialização das posições
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < Nx; i++) x[i] = i * dx;
    
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < Ny; j++) y[j] = j * dy;
    
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < Nz; k++) z[k] = k * dz;

    // Inicialização com zeros
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++)
                u[i*Ny*Nz + j*Nz + k] = 0.0;

    // Condição inicial gaussiana
    double cx = Lx/2, cy = Ly/2, cz = Lz/2;
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) {
                double dist2 = pow(x[i] - cx, 2) + pow(y[j] - cy, 2) + pow(z[k] - cz, 2);
                u[i*Ny*Nz + j*Nz + k] = exp(-dist2 / (2 * SIGMA * SIGMA));
            }

    // Cálculo das estatísticas iniciais
    compute_stats(u, &stats);
    print_stats(&stats, 0);

    struct timeval t0, t1;
    gettimeofday(&t0, NULL);

    for (int n = 0; n < nt; n++) {
        // Condições de contorno
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++) {
                u_new[i*Ny*Nz + j*Nz + 0] = 0.0;
                u_new[i*Ny*Nz + j*Nz + (Nz-1)] = 0.0;
            }
        
        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < Nx; i++)
            for (int k = 0; k < Nz; k++) {
                u_new[i*Ny*Nz + 0*Nz + k] = 0.0;
                u_new[i*Ny*Nz + (Ny-1)*Nz + k] = 0.0;
            }
        
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) {
                u_new[0*Ny*Nz + j*Nz + k] = 0.0;
                u_new[(Nx-1)*Ny*Nz + j*Nz + k] = 0.0;
            }

        // Evolução
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 1; i < Nx-1; i++)
            for (int j = 1; j < Ny-1; j++)
                for (int k = 1; k < Nz-1; k++)
                    u_new[i*Ny*Nz + j*Nz + k] = u[i*Ny*Nz + j*Nz + k] + NU * DT * (
                        (u[(i+1)*Ny*Nz + j*Nz + k] - 2*u[i*Ny*Nz + j*Nz + k] + u[(i-1)*Ny*Nz + j*Nz + k]) / (dx*dx) +
                        (u[i*Ny*Nz + (j+1)*Nz + k] - 2*u[i*Ny*Nz + j*Nz + k] + u[i*Ny*Nz + (j-1)*Nz + k]) / (dy*dy) +
                        (u[i*Ny*Nz + j*Nz + (k+1)] - 2*u[i*Ny*Nz + j*Nz + k] + u[i*Ny*Nz + j*Nz + (k-1)]) / (dz*dz)
                    );

        // Atualização
        #pragma omp parallel for collapse(3) schedule(static)
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++)
                for (int k = 0; k < Nz; k++)
                    u[i*Ny*Nz + j*Nz + k] = u_new[i*Ny*Nz + j*Nz + k];

        // Cálculo periódico de estatísticas
        if ((n+1) % (nt/10) == 0 || n == nt-1) {
            compute_stats(u, &stats);
            print_stats(&stats, n+1);
        }
    }

    gettimeofday(&t1, NULL);

    free(u);
    free(u_new);
    free(x);
    free(y);
    free(z);

    printf("Simulação 3D concluída (OpenMP).\n");
    printf("Tempo de execução: %.6lf segundos\n", elapsed(t0, t1));
    printf("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
    printf("Número máximo de threads: %d\n", omp_get_max_threads());
    
    return 0;
}