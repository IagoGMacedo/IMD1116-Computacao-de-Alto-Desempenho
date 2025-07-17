#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define IDX(i, j, N) ((i)*(N)+(j))

void exchange_borders(double *u, int N, int local_rows, int rank, int size, MPI_Request reqs[4]) {
    for(int i = 0; i < 4; i++) reqs[i] = MPI_REQUEST_NULL;
    
    if (rank > 0) { 
        MPI_Irecv(u, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &reqs[0]);
        MPI_Isend(u+N, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &reqs[1]);
    }
    if (rank < size-1) { 
        MPI_Irecv(u+(local_rows+1)*N, N, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &reqs[2]);
        MPI_Isend(u+local_rows*N, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &reqs[3]);
    }
}

void wait_borders(int rank, int size, MPI_Request reqs[4]) {
    if (rank > 0) {
        MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
        MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
    }
    if (rank < size-1) {
        MPI_Wait(&reqs[2], MPI_STATUS_IGNORE);
        MPI_Wait(&reqs[3], MPI_STATUS_IGNORE);
    }
}

void step(double *u, double *unew, int N, int local_rows, double alpha, double dt, double dx2) {
    #pragma omp parallel for
    for (int i = 1; i <= local_rows; i++) {
        for (int j = 1; j < N-1; j++) {
            unew[IDX(i,j,N)] = u[IDX(i,j,N)] +
                alpha * dt * (
                    (u[IDX(i-1,j,N)] + u[IDX(i+1,j,N)] +
                     u[IDX(i,j-1,N)] + u[IDX(i,j+1,N)] -
                     4*u[IDX(i,j,N)]) / dx2
                );
        }
    }
}

int main(int argc, char *argv[]) {
    int N = 1024;
    int STEPS = 1000;
    double alpha = 0.01, dt = 0.1, dx = 1.0, dx2 = dx*dx;

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc > 1) N = atoi(argv[1]);
    int local_rows = N / size;
    if (N % size != 0) {
        if (rank == 0) printf("N deve ser divisível por size\n");
        MPI_Finalize();
        return 1;
    }

    double max_dt = dx2 / (4.0 * alpha);
    if (rank == 0 && dt >= max_dt) {
        printf("Aviso: dt instável! Use dt < %.5f para estabilidade\n", max_dt);
    }

    double *u = calloc((local_rows+2)*N, sizeof(double));
    double *unew = calloc((local_rows+2)*N, sizeof(double));

    int global_center = N / 2;
    int center_rank = global_center / local_rows;
    if (rank == center_rank) {
        int local_center = global_center % local_rows;
        u[IDX(local_center + 1, global_center, N)] = 100.0;
    }

    double t0 = MPI_Wtime();
    for (int s = 0; s < STEPS; s++) {
        MPI_Request reqs[4];
        exchange_borders(u, N, local_rows, rank, size, reqs);

        #pragma omp parallel for
        for (int i = 2; i <= local_rows-1; i++) {
            for (int j = 1; j < N-1; j++) {
                unew[IDX(i,j,N)] = u[IDX(i,j,N)] +
                    alpha * dt * (
                        (u[IDX(i-1,j,N)] + u[IDX(i+1,j,N)] +
                         u[IDX(i,j-1,N)] + u[IDX(i,j+1,N)] -
                         4*u[IDX(i,j,N)]) / dx2
                    );
            }
        }

        wait_borders(rank, size, reqs);
        
        if (rank > 0) {
            #pragma omp parallel for
            for (int j = 1; j < N-1; j++) {
                unew[IDX(1,j,N)] = u[IDX(1,j,N)] +
                    alpha * dt * (
                        (u[IDX(0,j,N)] + u[IDX(2,j,N)] +
                         u[IDX(1,j-1,N)] + u[IDX(1,j+1,N)] -
                         4*u[IDX(1,j,N)]) / dx2
                    );
            }
        }
        
        if (rank < size-1) {
            #pragma omp parallel for
            for (int j = 1; j < N-1; j++) {
                unew[IDX(local_rows,j,N)] = u[IDX(local_rows,j,N)] +
                    alpha * dt * (
                        (u[IDX(local_rows-1,j,N)] + u[IDX(local_rows+1,j,N)] +
                         u[IDX(local_rows,j-1,N)] + u[IDX(local_rows,j+1,N)] -
                         4*u[IDX(local_rows,j,N)]) / dx2
                    );
            }
        }


        double *tmp = u;
        u = unew;
        unew = tmp;
    }
    double t1 = MPI_Wtime();

    if (rank == 0) {
        printf("N=%d, size=%d, tempo=%.3fs\n", N, size, t1-t0);
    }

    free(u);
    free(unew);
    MPI_Finalize();
    return 0;
}