#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
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
__global__ void set_borders_zero(double* u_new, int nx, int ny, int nz) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    // Face X = 0 e X = nx-1
    if (i < ny && j < nz) {
        u_new[0 * ny * nz + i * nz + j] = 0.0;
        u_new[(nx-1) * ny * nz + i * nz + j] = 0.0;
    }

    // Face Y = 0 e Y = ny-1
    if (i < nx && j < nz) {
        u_new[i * ny * nz + 0 * nz + j] = 0.0;
        u_new[i * ny * nz + (ny-1) * nz + j] = 0.0;
    }

    // Face Z = 0 e Z = nz-1
    if (i < nx && j < ny) {
        u_new[i * ny * nz + j * nz + 0] = 0.0;
        u_new[i * ny * nz + j * nz + (nz-1)] = 0.0;
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

// Kernel para redução de estatísticas usando memória compartilhada
__global__ void reduce_stats_kernel(double* u, double* stats, int size) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x * 4 + threadIdx.x;
    
    // Inicializa valores locais
    double my_sum = 0.0;
    double my_max = -DBL_MAX;
    double my_min = DBL_MAX;
    double my_sum_sq = 0.0;
    
    // Carrega até 4 elementos por thread
    if (i < size) {
        my_sum = u[i];
        my_max = u[i];
        my_min = u[i];
        my_sum_sq = u[i] * u[i];
    }
    if (i + blockDim.x < size) {
        double val = u[i + blockDim.x];
        my_sum += val;
        my_max = fmax(my_max, val);
        my_min = fmin(my_min, val);
        my_sum_sq += val * val;
    }
    if (i + 2*blockDim.x < size) {
        double val = u[i + 2*blockDim.x];
        my_sum += val;
        my_max = fmax(my_max, val);
        my_min = fmin(my_min, val);
        my_sum_sq += val * val;
    }
    if (i + 3*blockDim.x < size) {
        double val = u[i + 3*blockDim.x];
        my_sum += val;
        my_max = fmax(my_max, val);
        my_min = fmin(my_min, val);
        my_sum_sq += val * val;
    }
    
    // Memória compartilhada para 4 valores por thread (soma, max, min, sum_sq)
    int idx = tid * 4;
    sdata[idx] = my_sum;
    sdata[idx+1] = my_max;
    sdata[idx+2] = my_min;
    sdata[idx+3] = my_sum_sq;
    __syncthreads();
    
    // Redução em árvore
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            int sidx = tid * 4;
            int sidx2 = (tid + s) * 4;
            
            sdata[sidx] += sdata[sidx2];       // Soma
            sdata[sidx+1] = fmax(sdata[sidx+1], sdata[sidx2+1]); // Máximo
            sdata[sidx+2] = fmin(sdata[sidx+2], sdata[sidx2+2]); // Mínimo
            sdata[sidx+3] += sdata[sidx2+3];   // Soma dos quadrados
        }
        __syncthreads();
    }
    
    // Thread 0 escreve o resultado do bloco
    if (tid == 0) {
        int bidx = blockIdx.x * 4;
        stats[bidx] = sdata[0];       // Soma
        stats[bidx+1] = sdata[1];     // Máximo
        stats[bidx+2] = sdata[2];     // Mínimo
        stats[bidx+3] = sdata[3];     // Soma dos quadrados
    }
}

// Função wrapper para cálculo de estatísticas na GPU
void compute_stats_gpu(double* d_u, double* h_stats, int size) {
    int block_size = 256;
    int grid_size = (size + 4 * block_size - 1) / (4 * block_size);
    
    double* d_stats;
    CHECK(cudaMalloc(&d_stats, grid_size * 4 * sizeof(double)));
    
    size_t shared_mem_size = block_size * 4 * sizeof(double);
    reduce_stats_kernel<<<grid_size, block_size, shared_mem_size>>>(d_u, d_stats, size);
    CHECK(cudaGetLastError());
    
    double* block_stats = (double*)malloc(grid_size * 4 * sizeof(double));
    CHECK(cudaMemcpy(block_stats, d_stats, grid_size * 4 * sizeof(double), cudaMemcpyDeviceToHost));
    
    // Redução final na CPU (pequena)
    double total_sum = 0.0;
    double global_max = -DBL_MAX;
    double global_min = DBL_MAX;
    double total_sum_sq = 0.0;
    
    for (int i = 0; i < grid_size; i++) {
        int idx = i * 4;
        total_sum += block_stats[idx];
        global_max = fmax(global_max, block_stats[idx+1]);
        global_min = fmin(global_min, block_stats[idx+2]);
        total_sum_sq += block_stats[idx+3];
    }
    
    h_stats[0] = total_sum;
    h_stats[1] = global_max;
    h_stats[2] = global_min;
    h_stats[3] = total_sum_sq;
    
    free(block_stats);
    CHECK(cudaFree(d_stats));
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

    // Inicialização
    for (i = 0; i < Nx; i++) x[i] = i * dx;
    for (j = 0; j < Ny; j++) y[j] = j * dy;
    for (k = 0; k < Nz; k++) z[k] = k * dz;

    // Inicialização com zeros
    memset(u, 0, total_size * sizeof(double));
    memset(u_new, 0, total_size * sizeof(double));

    // Condição inicial gaussiana
    double cx = Lx/2, cy = Ly/2, cz = Lz/2;
    for (i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            for (k = 0; k < Nz; k++) {
                double dist2 = pow(x[i] - cx, 2) + pow(y[j] - cy, 2) + pow(z[k] - cz, 2);
                u[i*Ny*Nz + j*Nz + k] = exp(-dist2 / (2 * SIGMA * SIGMA));
            }
        }
    }

    // Alocação na GPU
    double *d_u, *d_u_new;
    CHECK(cudaMalloc(&d_u, total_size * sizeof(double)));
    CHECK(cudaMalloc(&d_u_new, total_size * sizeof(double)));
    
    // Copiar dados iniciais para GPU
    CHECK(cudaMemcpy(d_u, u, total_size * sizeof(double), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_u_new, u_new, total_size * sizeof(double), cudaMemcpyHostToDevice));

    // Configuração de kernels
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

    // Variáveis para estatísticas
    double stats[4];
    
    // Cálculo das estatísticas iniciais na GPU
    compute_stats_gpu(d_u, stats, total_size);
    print_stats(stats[0], stats[1], stats[2], stats[3], 0);

    struct timeval t0, t1;
    gettimeofday(&t0, NULL);

    // Loop de evolução temporal
    for (n = 0; n < nt; n++) {
        // Zerar bordas na GPU
        set_borders_zero<<<grid_borders, block_borders>>>(d_u_new, Nx, Ny, Nz);
        CHECK(cudaGetLastError());

        // Calcular evolução
        evolve_kernel<<<grid_evolve, block_evolve>>>(d_u, d_u_new, nu_dt, dx2, dy2, dz2, Nx, Ny, Nz);
        CHECK(cudaGetLastError());

        // Trocar ponteiros
        double* temp = d_u;
        d_u = d_u_new;
        d_u_new = temp;

        // Calcular estatísticas periodicamente na GPU
        if ((n+1) % (nt/10) == 0 || n == nt-1) {
            compute_stats_gpu(d_u, stats, total_size);
            print_stats(stats[0], stats[1], stats[2], stats[3], n+1);
        }
    }

    gettimeofday(&t1, NULL);

    // Limpeza
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