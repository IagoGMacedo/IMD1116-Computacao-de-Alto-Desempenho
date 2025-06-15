#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Parâmetros da simulação
#define N 10000000     // Número total de células (incluindo bordas)
#define NUM_STEPS 1000 // Número de passos de tempo
#define ALPHA 1.0      // Coeficiente de difusão
#define L 1.0          // Comprimento da barra

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Calcular divisão do domínio
    int total_inner = N - 2; // Células internas (excluindo bordas globais)
    int nlocal = total_inner / size;
    int rest = total_inner % size;
    int local_cells = nlocal + (rank < rest ? 1 : 0);

    // Alocar memória com células fantasmas
    double *u_old = (double *)malloc((local_cells + 2) * sizeof(double));
    double *u_new = (double *)malloc((local_cells + 2) * sizeof(double));

    // Inicializar vetores locais
    for (int i = 0; i < local_cells + 2; i++)
    {
        u_old[i] = 0.0;
        u_new[i] = 0.0;
    }

    // Configurar condições de contorno globais
    if (rank == 0)
    {
        u_old[0] = 100.0; // Borda esquerda global
    }
    if (rank == size - 1)
    {
        u_old[local_cells + 1] = 0.0; // Borda direita global
    }

    // Tamanho do passo espacial e temporal
    double dx = L / (N - 1);
    double dt = 0.1 * (dx * dx) / ALPHA;

    // Identificar vizinhos
    int left_rank = (rank > 0) ? rank - 1 : MPI_PROC_NULL;
    int right_rank = (rank < size - 1) ? rank + 1 : MPI_PROC_NULL;

    double start_time = MPI_Wtime();

    for (int step = 0; step < NUM_STEPS; step++)
    {
        MPI_Request requests[4];
        int request_count = 0;

        // Postar recebimentos primeiro para evitar deadlocks
        if (left_rank != MPI_PROC_NULL)
        {
            MPI_Irecv(&u_old[0], 1, MPI_DOUBLE, left_rank, 0,
                      MPI_COMM_WORLD, &requests[request_count++]);
        }
        if (right_rank != MPI_PROC_NULL)
        {
            MPI_Irecv(&u_old[local_cells + 1], 1, MPI_DOUBLE, right_rank, 1,
                      MPI_COMM_WORLD, &requests[request_count++]);
        }

        // Postar envios
        if (right_rank != MPI_PROC_NULL)
        {
            MPI_Isend(&u_old[local_cells], 1, MPI_DOUBLE, right_rank, 0,
                      MPI_COMM_WORLD, &requests[request_count++]);
        }
        if (left_rank != MPI_PROC_NULL)
        {
            MPI_Isend(&u_old[1], 1, MPI_DOUBLE, left_rank, 1,
                      MPI_COMM_WORLD, &requests[request_count++]);
        }

        // Esperar por todas as operações de comunicação
        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);

        // Atualizar todas as células locais
        for (int i = 1; i <= local_cells; i++)
        {
            u_new[i] = u_old[i] + ALPHA * dt / (dx * dx) *
                                      (u_old[i - 1] - 2 * u_old[i] + u_old[i + 1]);
        }

        // Trocar buffers para próxima iteração
        double *temp = u_old;
        u_old = u_new;
        u_new = temp;
    }

    double local_time = MPI_Wtime() - start_time;
    double max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("Tempo da versão 2 (Async/Wait): %.6f segundos\n", max_time);
    }

    free(u_old);
    free(u_new);
    MPI_Finalize();
    return 0;
}