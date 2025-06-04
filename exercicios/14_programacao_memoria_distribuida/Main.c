//mpicc -O2 -o programa Main.c
//mpirun -np 2 ./programa


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;
    int n_iter = 10000;
    int sizes[] = {8, 64, 512, 4096, 32768, 262144, 1048576}; // 8B, 64B, 512B, 4KB, 32KB, 256KB, 1MB
    int n_sizes = sizeof(sizes) / sizeof(sizes[0]);
    char *buffer;
    double t_start, t_end, t_total;
    int i, s;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 2) {
        if (rank == 0) printf("Este programa deve ser executado com exatamente 2 processos.\n");
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) printf("Tamanho(bytes),Tempo_total(s),Tempo_medio(s)\n");

    for (s = 0; s < n_sizes; s++) {
        int msg_size = sizes[s];
        buffer = (char*) malloc(msg_size);

        memset(buffer, 1, msg_size);

        MPI_Barrier(MPI_COMM_WORLD);

        t_start = MPI_Wtime();

        for (i = 0; i < n_iter; i++) {
            if (rank == 0) {
                MPI_Send(buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
                MPI_Recv(buffer, msg_size, MPI_CHAR, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else if (rank == 1) {
                MPI_Recv(buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(buffer, msg_size, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
            }
        }

        t_end = MPI_Wtime();
        t_total = t_end - t_start;

        if (rank == 0) {
            double t_medio = t_total / (2 * n_iter);
            printf("%d,%.6f,%.9f\n", msg_size, t_total, t_medio);
        }

        free(buffer);
    }

    MPI_Finalize();
    return 0;
}