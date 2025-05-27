#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define Lx 1.0
#define Ly 1.0
#define Lz 1.0
#define Nx 181
#define Ny 181
#define Nz 181
#define NU 0.01
#define DT 0.0005
#define T 0.1
#define SIGMA 0.1

double elapsed(struct timeval t0, struct timeval t1) {
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1e6;
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
    int i, j, k, n;

    for (i = 0; i < Nx; i++) x[i] = i * dx;
    for (j = 0; j < Ny; j++) y[j] = j * dy;
    for (k = 0; k < Nz; k++) z[k] = k * dz;

    for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
            for (k = 0; k < Nz; k++)
                u[i*Ny*Nz + j*Nz + k] = 0.0;

    FILE *fp = fopen("evolucao3d.csv", "w");
    if (!fp) {
        printf("Erro ao abrir arquivo!\n");
        return 1;
    }
    for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
            for (k = 0; k < Nz; k++)
                fprintf(fp, "%lf%c", u[i*Ny*Nz + j*Nz + k], (k < Nz-1) ? ',' : '\n');
    fprintf(fp, "\n");

    double cx = Lx/2, cy = Ly/2, cz = Lz/2;
    for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
            for (k = 0; k < Nz; k++) {
                double dist2 = pow(x[i] - cx, 2) + pow(y[j] - cy, 2) + pow(z[k] - cz, 2);
                u[i*Ny*Nz + j*Nz + k] = exp(-dist2 / (2 * SIGMA * SIGMA));
            }

    for (i = 0; i < Nx; i++)
        for (j = 0; j < Ny; j++)
            for (k = 0; k < Nz; k++)
                fprintf(fp, "%lf%c", u[i*Ny*Nz + j*Nz + k], (k < Nz-1) ? ',' : '\n');
    fprintf(fp, "\n");

    struct timeval t0, t1;
    gettimeofday(&t0, NULL);

    for (n = 0; n < nt; n++) {
        for (i = 0; i < Nx; i++)
            for (j = 0; j < Ny; j++) {
                u_new[i*Ny*Nz + j*Nz + 0] = 0.0;
                u_new[i*Ny*Nz + j*Nz + (Nz-1)] = 0.0;
            }
        for (i = 0; i < Nx; i++)
            for (k = 0; k < Nz; k++) {
                u_new[i*Ny*Nz + 0*Nz + k] = 0.0;
                u_new[i*Ny*Nz + (Ny-1)*Nz + k] = 0.0;
            }
        for (j = 0; j < Ny; j++)
            for (k = 0; k < Nz; k++) {
                u_new[0*Ny*Nz + j*Nz + k] = 0.0;
                u_new[(Nx-1)*Ny*Nz + j*Nz + k] = 0.0;
            }

        for (i = 1; i < Nx-1; i++)
            for (j = 1; j < Ny-1; j++)
                for (k = 1; k < Nz-1; k++)
                    u_new[i*Ny*Nz + j*Nz + k] = u[i*Ny*Nz + j*Nz + k] + NU * DT * (
                        (u[(i+1)*Ny*Nz + j*Nz + k] - 2*u[i*Ny*Nz + j*Nz + k] + u[(i-1)*Ny*Nz + j*Nz + k]) / (dx*dx) +
                        (u[i*Ny*Nz + (j+1)*Nz + k] - 2*u[i*Ny*Nz + j*Nz + k] + u[i*Ny*Nz + (j-1)*Nz + k]) / (dy*dy) +
                        (u[i*Ny*Nz + j*Nz + (k+1)] - 2*u[i*Ny*Nz + j*Nz + k] + u[i*Ny*Nz + j*Nz + (k-1)]) / (dz*dz)
                    );

        for (i = 0; i < Nx; i++)
            for (j = 0; j < Ny; j++)
                for (k = 0; k < Nz; k++)
                    u[i*Ny*Nz + j*Nz + k] = u_new[i*Ny*Nz + j*Nz + k];

        if ((n+1) % (nt/4) == 0 || n == nt-1) {
            for (i = 0; i < Nx; i++)
                for (j = 0; j < Ny; j++)
                    for (k = 0; k < Nz; k++)
                        fprintf(fp, "%lf%c", u[i*Ny*Nz + j*Nz + k], (k < Nz-1) ? ',' : '\n');
            fprintf(fp, "\n");
        }
    }

    gettimeofday(&t1, NULL);

    fclose(fp);
    free(u);
    free(u_new);
    free(x);
    free(y);
    free(z);

    printf("Simulação 3D concluída (codigoSequencial). Resultados salvos em evolucao3d.csv\n");
    printf("Tempo de execução: %.6lf segundos\n", elapsed(t0, t1));
    printf("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
    return 0;
}