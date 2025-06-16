#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define ROWS_ARG_IDX 1
#define COLS_ARG_IDX 2

int main(int argc, char *argv[]) {
    int processRank, numberOfProcesses;
    int rowNumber, colNumber;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    if (processRank == 0) {
        rowNumber = atoi(argv[ROWS_ARG_IDX]);
        colNumber = atoi(argv[COLS_ARG_IDX]);
    }

    MPI_Bcast(&rowNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&colNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);


    int colsPerProcess = colNumber / numberOfProcesses;
    double *matrix = NULL;
    double *vectorX = NULL;
    double *resultVector = NULL;

    MPI_Datatype block_of_columns_type, resized_block_type;
    MPI_Type_vector(rowNumber, colsPerProcess, colNumber, MPI_DOUBLE, &block_of_columns_type);
    MPI_Type_create_resized(block_of_columns_type, 0, colsPerProcess * sizeof(double), &resized_block_type);
    MPI_Type_commit(&resized_block_type);

    double *localBlock = (double *)malloc(rowNumber * colsPerProcess * sizeof(double));
    double *localX = (double *)malloc(colsPerProcess * sizeof(double));
    double *y_local = (double *)calloc(rowNumber, sizeof(double)); 

    if (processRank == 0) {
        matrix = (double *)malloc(rowNumber * colNumber * sizeof(double));
        vectorX = (double *)malloc(colNumber * sizeof(double));
        resultVector = (double *)malloc(rowNumber * sizeof(double));

        srand(time(NULL));
        for (int i = 0; i < rowNumber * colNumber; i++)
            matrix[i] = rand() % 10;

        for (int i = 0; i < colNumber; i++)
            vectorX[i] = rand() % 10;
    }

    double start_time, end_time, elapsed_time;
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    MPI_Scatter(matrix, 1, resized_block_type, localBlock, rowNumber * colsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Scatter(vectorX, colsPerProcess, MPI_DOUBLE, localX, colsPerProcess, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < rowNumber; i++) {
        for (int j = 0; j < colsPerProcess; j++) {
            y_local[i] += localBlock[i * colsPerProcess + j] * localX[j];
        }
    }

    MPI_Reduce(y_local, resultVector, rowNumber, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    end_time = MPI_Wtime();
    elapsed_time = end_time - start_time;

    MPI_Type_free(&block_of_columns_type);
    MPI_Type_free(&resized_block_type);

    if (processRank == 0) {
        printf("Tempo decorrido: %f segundos\n", elapsed_time);
        printf("Cálculo executado com %d linhas, %d colunas e %d processos.\n",
               rowNumber, colNumber, numberOfProcesses);

        free(matrix);
        free(vectorX);
        free(resultVector);
    }

    free(localBlock);
    free(localX);
    free(y_local);

    MPI_Finalize();
    return 0;
}