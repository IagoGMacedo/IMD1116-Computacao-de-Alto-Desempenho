#include <stdio.h>
#include <string.h>
#include <cuda_runtime.h>

#define BLOCK_SIZE 8
#define HALO 1

__global__ void atualiza(double *vnew, double *vold, int nx, int ny, int nz, double alpha)
{
    // Coordenadas do bloco sem halo
    int bx = blockIdx.x * BLOCK_SIZE;
    int by = blockIdx.y * BLOCK_SIZE;
    int bz = blockIdx.z * BLOCK_SIZE;
    
    // Coordenadas locais da thread
    int tx = threadIdx.x;
    int ty = threadIdx.y;
    int tz = threadIdx.z;
    
    // Memória compartilhada com halo (BLOCK_SIZE + 2*HALO)
    __shared__ double s_data[BLOCK_SIZE+2][BLOCK_SIZE+2][BLOCK_SIZE+2];
    
    // Coordenadas globais para carregamento
    int gx = bx + tx - HALO;
    int gy = by + ty - HALO;
    int gz = bz + tz - HALO;
    
    // Carregar bloco central + halo na memória compartilhada
    if (gx >= 0 && gx < nx && gy >= 0 && gy < ny && gz >= 0 && gz < nz) {
        s_data[tz][ty][tx] = vold[gz * ny * nx + gy * nx + gx];
    } else {
        s_data[tz][ty][tx] = 0.0; // Condições de contorno implícitas
    }
    
    __syncthreads();
    
    // Apenas threads centrais calculam (não fazem parte do halo)
    if (tx >= HALO && tx < BLOCK_SIZE+HALO && 
        ty >= HALO && ty < BLOCK_SIZE+HALO && 
        tz >= HALO && tz < BLOCK_SIZE+HALO) {
        
        // Coordenadas globais para escrita
        gx = bx + tx - HALO;
        gy = by + ty - HALO;
        gz = bz + tz - HALO;
        
        if (gx > 0 && gx < nx-1 && gy > 0 && gy < ny-1 && gz > 0 && gz < nz-1) {
            double val = s_data[tz][ty][tx];
            double sum = s_data[tz][ty][tx+1] + s_data[tz][ty][tx-1] +
                         s_data[tz][ty+1][tx] + s_data[tz][ty-1][tx] +
                         s_data[tz+1][ty][tx] + s_data[tz-1][ty][tx];
            
            vnew[gz * ny * nx + gy * nx + gx] = val + alpha * (sum - 6.0 * val);
        }
    }
}

int main()
{
    // Configuração do problema
    const int nx = 381, ny = 381, nz = 381;
    const int nt = 381;
    const double alpha = 0.1;
    
    // Alocar e inicializar memória na CPU
    double *h_vold = (double*)malloc(nx * ny * nz * sizeof(double));
    double *h_result = (double*)malloc(nx * ny * nz * sizeof(double));
    memset(h_vold, 0, nx * ny * nz * sizeof(double));
    h_vold[(nz/2)*ny*nx + (ny/2)*nx + (nx/2)] = 1.0; // Fonte no centro
    
    // Alocar memória na GPU
    double *d_vold, *d_vnew;
    size_t size = nx * ny * nz * sizeof(double);
    cudaMalloc(&d_vold, size);
    cudaMalloc(&d_vnew, size);
    cudaMemcpy(d_vold, h_vold, size, cudaMemcpyHostToDevice);
    
    // Configurar kernel - blocos de 10x10x10 threads (8x8x8 úteis + halo)
    dim3 threads(BLOCK_SIZE+2, BLOCK_SIZE+2, BLOCK_SIZE+2);
    dim3 grid(
        (nx + BLOCK_SIZE - 1) / BLOCK_SIZE,
        (ny + BLOCK_SIZE - 1) / BLOCK_SIZE,
        (nz + BLOCK_SIZE - 1) / BLOCK_SIZE
    );
    
    // Medir tempo
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    
    // Executar iterações
    for (int t = 0; t < nt; t++) {
        atualiza<<<grid, threads>>>(d_vnew, d_vold, nx, ny, nz, alpha);
        std::swap(d_vnew, d_vold);
    }
    
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms;
    cudaEventElapsedTime(&ms, start, stop);
    printf("Tempo de execução: %.2f ms\n", ms);
    
    // Copiar resultados e limpar
    cudaMemcpy(h_result, d_vold, size, cudaMemcpyDeviceToHost);
    
    // Salvar resultados (opcional)
    FILE *file = fopen("resultados_otimizados.txt", "w");
    if (file) {
        for (int i = 0; i < nx*ny*nz; i++) {
            fprintf(file, "%.6e\n", h_result[i]);
        }
        fclose(file);
    }
    
    // Liberar recursos
    cudaFree(d_vold);
    cudaFree(d_vnew);
    free(h_vold);
    free(h_result);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    
    return 0;
}