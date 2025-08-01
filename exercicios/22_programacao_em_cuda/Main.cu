#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

#define Lx 1.0
#define Ly 1.0
#define Lz 1.0
#define Nx 301
#define Ny 301
#define Nz 301
#define NU 0.1
#define DT 0.00002
#define T 0.01
#define SIGMA 0.1

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

__global__ void set_borders_zero(double* u_new, int nx, int ny, int nz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < ny && j < nz) {
        u_new[0 * ny * nz + i * nz + j] = 0.0;
        u_new[(nx-1) * ny * nz + i * nz + j] = 0.0;
    }

    if (i < nx && j < nz) {
        u_new[i * ny * nz + 0 * nz + j] = 0.0;
        u_new[i * ny * nz + (ny-1) * nz + j] = 0.0;
    }

    if (i < nx && j < ny) {
        u_new[i * ny * nz + j * nz + 0] = 0.0;
        u_new[i * ny * nz + j * nz + (nz-1)] = 0.0;
    }
}

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

// Kernel com redução simplificada (sem shared memory)
__global__ void partial_stats_kernel(double* u, double* partial_sums, double* partial_maxs, 
                                    double* partial_mins, double* partial_sum_sqs, int size, int partial_size) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    double local_sum = 0.0;
    double local_max = -DBL_MAX;
    double local_min = DBL_MAX;
    double local_sum_sq = 0.0;

    for (int i = idx; i < size; i += stride) {
        double val = u[i];
        local_sum += val;
        local_max = fmax(local_max, val);
        local_min = fmin(local_min, val);
        local_sum_sq += val * val;
    }

    // Apenas threads dentro do limite de partial_size escrevem
    if (idx < partial_size) {
        partial_sums[idx] = local_sum;
        partial_maxs[idx] = local_max;
        partial_mins[idx] = local_min;
        partial_sum_sqs[idx] = local_sum_sq;
    }
}

void compute_stats_gpu(double* d_u, double* h_stats, int size) {
    int block_size = 256;
    int grid_size = (size + block_size - 1) / block_size;
    if (grid_size > 1024) grid_size = 1024;

    // Tamanho real do array parcial
    int partial_size = grid_size * block_size;

    double *d_partial_sums, *d_partial_maxs, *d_partial_mins, *d_partial_sum_sqs;
    CHECK(cudaMalloc(&d_partial_sums, partial_size * sizeof(double)));
    CHECK(cudaMalloc(&d_partial_maxs, partial_size * sizeof(double)));
    CHECK(cudaMalloc(&d_partial_mins, partial_size * sizeof(double)));
    CHECK(cudaMalloc(&d_partial_sum_sqs, partial_size * sizeof(double)));

    partial_stats_kernel<<<grid_size, block_size>>>(d_u, d_partial_sums, d_partial_maxs, 
                                                  d_partial_mins, d_partial_sum_sqs, size, partial_size);
    CHECK(cudaDeviceSynchronize());

    double* partial_sums = (double*)malloc(partial_size * sizeof(double));
    double* partial_maxs = (double*)malloc(partial_size * sizeof(double));
    double* partial_mins = (double*)malloc(partial_size * sizeof(double));
    double* partial_sum_sqs = (double*)malloc(partial_size * sizeof(double));
    
    CHECK(cudaMemcpy(partial_sums, d_partial_sums, partial_size * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(partial_maxs, d_partial_maxs, partial_size * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(partial_mins, d_partial_mins, partial_size * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(partial_sum_sqs, d_partial_sum_sqs, partial_size * sizeof(double), cudaMemcpyDeviceToHost));

    double total_sum = 0.0;
    double global_max = -DBL_MAX;
    double global_min = DBL_MAX;
    double total_sum_sq = 0.0;

    for (int i = 0; i < partial_size; i++) {
        total_sum += partial_sums[i];
        global_max = fmax(global_max, partial_maxs[i]);
        global_min = fmin(global_min, partial_mins[i]);
        total_sum_sq += partial_sum_sqs[i];
    }

    h_stats[0] = total_sum;
    h_stats[1] = global_max;
    h_stats[2] = global_min;
    h_stats[3] = total_sum_sq;

    free(partial_sums);
    free(partial_maxs);
    free(partial_mins);
    free(partial_sum_sqs);
    CHECK(cudaFree(d_partial_sums));
    CHECK(cudaFree(d_partial_maxs));
    CHECK(cudaFree(d_partial_mins));
    CHECK(cudaFree(d_partial_sum_sqs));
}

void print_stats(double mass, double max_val, double min_val, double l2_squared, int step) {
    printf("Passo %d:\n", step);
    printf("  Massa total: %.6f\n", mass);
    printf("  Valor máximo: %.6f\n", max_val);
    printf("  Valor mínimo: %.6f\n", min_val);
    printf("  Norma L2: %.6f\n\n", sqrt(l2_squared));
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
    int total_size = Nx * Ny * Nz;

    double *u = (double*)malloc(total_size * sizeof(double));
    double *u_new = (double*)malloc(total_size * sizeof(double));
    double *x = (double*)malloc(Nx * sizeof(double));
    double *y = (double*)malloc(Ny * sizeof(double));
    double *z = (double*)malloc(Nz * sizeof(double));
    int i, j, k, n;

    for (i = 0; i < Nx; i++) x[i] = i * dx;
    for (j = 0; j < Ny; j++) y[j] = j * dy;
    for (k = 0; k < Nz; k++) z[k] = k * dz;

    memset(u, 0, total_size * sizeof(double));
    memset(u_new, 0, total_size * sizeof(double));

    double cx = Lx/2, cy = Ly/2, cz = Lz/2;
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                double dist2 = pow(x[i] - cx, 2) + pow(y[j] - cy, 2) + pow(z[k] - cz, 2);
                u[i*Ny*Nz + j*Nz + k] = exp(-dist2 / (2 * SIGMA * SIGMA));
            }
        }
    }

    double *d_u, *d_u_new;
    CHECK(cudaMalloc(&d_u, total_size * sizeof(double)));
    CHECK(cudaMalloc(&d_u_new, total_size * sizeof(double)));
    
    CHECK(cudaMemcpy(d_u, u, total_size * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_u_new, u_new, total_size * sizeof(double), cudaMemcpyHostToDevice));

    dim3 block_borders(16, 16);
    int max_dim1 = (Nx > Ny) ? Nx : Ny;
    int max_dim2 = (Ny > Nz) ? Ny : Nz;
    if (Nz > max_dim2) max_dim2 = Nz;
    dim3 grid_borders(
        (max_dim1 + block_borders.x - 1) / block_borders.x,
        (max_dim2 + block_borders.y - 1) / block_borders.y
    );
    
    dim3 block_evolve(8, 8, 4);
    dim3 grid_evolve(
        (Nx - 2 + block_evolve.x - 1) / block_evolve.x,
        (Ny - 2 + block_evolve.y - 1) / block_evolve.y,
        (Nz - 2 + block_evolve.z - 1) / block_evolve.z
    );

    double stats[4];
    compute_stats_gpu(d_u, stats, total_size);
    print_stats(stats[0], stats[1], stats[2], stats[3], 0);

    struct timeval t0, t1;
    gettimeofday(&t0, NULL);

    for (n = 0; n < nt; n++) {
        set_borders_zero<<<grid_borders, block_borders>>>(d_u_new, Nx, Ny, Nz);
        CHECK(cudaGetLastError());

        evolve_kernel<<<grid_evolve, block_evolve>>>(d_u, d_u_new, nu_dt, dx2, dy2, dz2, Nx, Ny, Nz);
        CHECK(cudaGetLastError());

        double* temp = d_u;
        d_u = d_u_new;
        d_u_new = temp;

        if ((n+1) % (nt/10) == 0 || n == nt-1) {
            compute_stats_gpu(d_u, stats, total_size);
            print_stats(stats[0], stats[1], stats[2], stats[3], n+1);
        }
    }

    gettimeofday(&t1, NULL);

    free(u);
    free(u_new);
    free(x);
    free(y);
    free(z);
    CHECK(cudaFree(d_u));
    CHECK(cudaFree(d_u_new));

    printf("Simulação 3D concluída (CUDA)\n");
    printf("Tempo de execução: %.6lf segundos\n", elapsed(t0, t1));
    printf("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
    
    return 0;
}