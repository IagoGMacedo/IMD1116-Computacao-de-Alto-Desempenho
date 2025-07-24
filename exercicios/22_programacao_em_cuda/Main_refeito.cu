#include <stdio.h>

__global__ void atualiza(double *vnew, double *vold, int nx, int ny, int nz, double alpha)
{
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    int idx = z * ny * nx + y * nx + x;
    if (x > 0 && x < nx - 1 && y > 0 && y < ny - 1 && z > 0 && z < nz - 1)
    {
        int xm = idx - 1;
        int xp = idx + 1;
        int ym = idx - nx;
        int yp = idx + nx;
        int zm = idx - nx * ny;
        int zp = idx + nx * ny;
        vnew[idx] = vold[idx] + alpha * (vold[xp] + vold[xm] +
                                         vold[yp] + vold[ym] +
                                         vold[zp] + vold[zm] - 6 * vold[idx]);
    }
}

int main()
{
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // declaração de variáveis
    int nx = 381, ny = 381, nz = 381;
    int nt = 381;
    double alpha = 0.1;

    // Aloca memória na CPU para os dados iniciais (h_vold) e resultados (h_result)
    double *h_vold = (double *)malloc(nx * ny * nz * sizeof(double));
    double *h_result = (double *)malloc(nx * ny * nz * sizeof(double));

    // Inicializa h_vold (exemplo: tudo zero, com um ponto central ativo)
    memset(h_vold, 0, nx * ny * nz * sizeof(double));
    h_vold[(nz / 2) * ny * nx + (ny / 2) * nx + (nx / 2)] = 1.0; // fonte no centro

    // Aloca memória na GPU
    double *d_vold, *d_vnew;
    int size = nx * ny * nz * sizeof(double);

    cudaMalloc(&d_vold, size);
    cudaMalloc(&d_vnew, size);

    // Copia dados iniciais para a GPU
    cudaMemcpy(d_vold, h_vold, size, cudaMemcpyHostToDevice);

    // Define bloco 3D
    int bx = 8, by = 8, bz = 8;
    dim3 threads(bx, by, bz);
    dim3 grid((nx + bx - 1) / bx, (ny + by - 1) / by, (nz + bz - 1) / bz);
    for (int t = 0; t < nt; t++)
    {
        atualiza<<<grid, threads>>>(d_vnew, d_vold, nx, ny, nz, alpha);
        double *tmp = d_vold;
        d_vold = d_vnew;
        d_vnew = tmp;
    }
    cudaMemcpy(h_result, d_vold, size, cudaMemcpyDeviceToHost);

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

    float ms;
    cudaEventElapsedTime(&ms, start, stop);
    printf("Tempo: %f ms\n", ms);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    FILE *file = fopen("resultados_gpu.txt", "w");
    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo!\n");
        return 1;
    }

    for (int i = 0; i < nx * ny * nz; i++)
    {
        fprintf(file, "%e\n", h_result[i]);
    }

    fclose(file); // Fecha o arquivo
    printf("Resultados salvos em 'resultados_gpu.txt'.\n");

    return 0;
}