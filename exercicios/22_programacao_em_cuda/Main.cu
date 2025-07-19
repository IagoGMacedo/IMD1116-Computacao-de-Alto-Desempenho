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

// Macro para verificação de erros CUDA
#define CHECK(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "Erro CUDA (arquivo: %s, linha: %d): %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

double elapsed(struct timeval t0, struct timeval t1) {
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec) / 1e6;
}

// Kernel para zerar as bordas
__global__ void set_borders_zero(double* u, int nx, int ny, int nz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    // Face x = 0
    if (i < ny && j < nz) {
        u[0 * ny * nz + i * nz + j] = 0.0;
    }
    
    // Face x = nx-1
    if (i < ny && j < nz) {
        u[(nx-1) * ny * nz + i * nz + j] = 0.0;
    }
    
    // Face y = 0
    if (i < nx && j < nz) {
        u[i * ny * nz + 0 * nz + j] = 0.0;
    }
    
    // Face y = ny-1
    if (i < nx && j < nz) {
        u[i * ny * nz + (ny-1) * nz + j] = 0.0;
    }
    
    // Face z = 0
    if (i < nx && j < ny) {
        u[i * ny * nz + j * nz + 0] = 0.0;
    }
    
    // Face z = nz-1
    if (i < nx && j < ny) {
        u[i * ny * nz + j * nz + (nz-1)] = 0.0;
    }
}

// Kernel principal de evolução
__global__ void evolve_kernel(double* u, double* u_new, 
                              double nu_dt, double dx2, double dy2, double dz2,
                              int nx, int ny, int nz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;
    int k = blockIdx.z * blockDim.z + threadIdx.z + 1;

    if (i < nx - 1 && j < ny - 1 && k < nz - 1) {
        int idx = i * ny * nz + j * nz + k;
        
        double u_ijk = u[idx];
        double d2x = (u[(i+1)*ny*nz + j*nz + k] - 2*u_ijk + u[(i-1)*ny*nz + j*nz + k]) / dx2;
        double d2y = (u[i*ny*nz + (j+1)*nz + k] - 2*u_ijk + u[i*ny*nz + (j-1)*nz + k]) / dy2;
        double d2z = (u[i*ny*nz + j*nz + (k+1)] - 2*u_ijk + u[i*ny*nz + j*nz + (k-1)]) / dz2;
        
        u_new[idx] = u_ijk + nu_dt * (d2x + d2y + d2z);
    }
}

int main() {
    int nt = (int)(T / DT);
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);
    double dz = Lz / (Nz - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double dz2 = dz * dz;
    double nu_dt = NU * DT;

    double *u = (double*)malloc(Nx * Ny * Nz * sizeof(double));
    double *u_new = (double*)malloc(Nx * Ny * Nz * sizeof(double));
    double *x = (double*)malloc(Nx * sizeof(double));
    double *y = (double*)malloc(Ny * sizeof(double));
    double *z = (double*)malloc(Nz * sizeof(double));
    int i, j, k, n;

    for (i = 0; i < Nx; i++) x[i] = i * dx;
    for (j = 0; j < Ny; j++) y[j] = j * dy;
    for (k = 0; k < Nz; k++) z[k] = k * dz;

    // Inicialização
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

    // Condição inicial gaussiana
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

    // Alocação na GPU
    double *d_u, *d_u_new;
    CHECK(cudaMalloc(&d_u, Nx * Ny * Nz * sizeof(double)));
    CHECK(cudaMalloc(&d_u_new, Nx * Ny * Nz * sizeof(double)));
    
    // Copiar dados iniciais para GPU
    CHECK(cudaMemcpy(d_u, u, Nx * Ny * Nz * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_u_new, u_new, Nx * Ny * Nz * sizeof(double), cudaMemcpyHostToDevice));

    // Configuração de kernels
    dim3 block_borders(16, 16);
    dim3 grid_borders((Ny + 15) / 16, (Nz + 15) / 16);
    
    dim3 block_evolve(8, 8, 4);
    dim3 grid_evolve(
        (Nx - 2 + block_evolve.x - 1) / block_evolve.x,
        (Ny - 2 + block_evolve.y - 1) / block_evolve.y,
        (Nz - 2 + block_evolve.z - 1) / block_evolve.z
    );

    struct timeval t0, t1;
    gettimeofday(&t0, NULL);

    for (n = 0; n < nt; n++) {
        // Zerar bordas na GPU
        set_borders_zero<<<grid_borders, block_borders>>>(d_u_new, Nx, Ny, Nz);
        CHECK(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());

        // Calcular evolução
        evolve_kernel<<<grid_evolve, block_evolve>>>(d_u, d_u_new, nu_dt, dx2, dy2, dz2, Nx, Ny, Nz);
        CHECK(cudaGetLastError());
        CHECK(cudaDeviceSynchronize());

        // Trocar ponteiros
        double* temp = d_u;
        d_u = d_u_new;
        d_u_new = temp;

        // Escrever resultados periodicamente
        if ((n+1) % (nt/4) == 0 || n == nt-1) {
            CHECK(cudaMemcpy(u, d_u, Nx * Ny * Nz * sizeof(double), cudaMemcpyDeviceToHost));
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
    CHECK(cudaFree(d_u));
    CHECK(cudaFree(d_u_new));

    printf("Simulação 3D concluída (CUDA). Resultados salvos em evolucao3d.csv\n");
    printf("Tempo de execução: %.6lf segundos\n", elapsed(t0, t1));
    printf("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
    
    return 0;
}