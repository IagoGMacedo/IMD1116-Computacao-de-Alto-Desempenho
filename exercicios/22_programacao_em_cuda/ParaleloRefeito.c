#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

void atualiza(double *vnew, double *vold, int nx, int ny, int nz, double alpha)
{
#pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
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
        }
    }
}

int main()
{
    double start = omp_get_wtime();

    // declaração de variáveis
    int nx = 381, ny = 381, nz = 381;
    int nt = 381;
    double alpha = 0.1;

    // Aloca memória para os dados iniciais (vold)
    double *vold = (double *)malloc(nx * ny * nz * sizeof(double));
    double *vnew = (double *)malloc(nx * ny * nz * sizeof(double));

#pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                int idx = z * ny * nx + y * nx + x;
                vold[idx] = 0.0;
                if (x == nx / 2 && y == ny / 2 && z == nz / 2)
                {
                    vold[idx] = 1.0;
                }
            }
        }
    }

    // Executa as iterações
    for (int t = 0; t < nt; t++)
    {
        atualiza(vnew, vold, nx, ny, nz, alpha);
        double *tmp = vold;
        vold = vnew;
        vnew = tmp;
    }

    double end = omp_get_wtime();
    printf("Tempo: %f ms\n", (end - start) * 1000);

    FILE *file = fopen("resultados_cpu.txt", "w");
    if (file == NULL)
    {
        printf("Erro ao abrir o arquivo!\n");
        return 1;
    }

    for (int i = 0; i < nx * ny * nz; i++)
    {
        fprintf(file, "%e\n", vold[i]);
    }

    fclose(file); // Fecha o arquivo
    printf("Resultados salvos em 'resultados_cpu.txt'.\n");

    free(vold);
    free(vnew);

    return 0;
}