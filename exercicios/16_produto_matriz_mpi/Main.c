    #include <stdio.h>
    #include <stdlib.h>
    #include <time.h>
    #include <mpi.h>

    #define ROWS_ARG_IDX 1
    #define COLS_ARG_IDX 2

    int main(int argc, char *argv[])
    {
        int processRank, numberOfProcesses;
        int rowNumber, colNumber;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &processRank);
        MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

        if (processRank == 0)
        {
            rowNumber = atoi(argv[ROWS_ARG_IDX]);
            colNumber = atoi(argv[COLS_ARG_IDX]);
        }

        MPI_Bcast(&rowNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&colNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);

        double *matrix = NULL;
        double *vectorX = (double *)malloc(colNumber * sizeof(double));
        double *resultVector = NULL;

        int rowsPerProcess = rowNumber / numberOfProcesses;
        double *localMatrix = (double *)malloc(rowsPerProcess * colNumber * sizeof(double));
        double *localResult = (double *)malloc(rowsPerProcess * sizeof(double));

        if (processRank == 0)
        {
            matrix = (double *)malloc(rowNumber * colNumber * sizeof(double));
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

        MPI_Scatter(matrix, rowsPerProcess * colNumber, MPI_DOUBLE,
                    localMatrix, rowsPerProcess * colNumber, MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        MPI_Bcast(vectorX, colNumber, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int rowIdx = 0; rowIdx < rowsPerProcess; rowIdx++)
        {
            localResult[rowIdx] = 0.0;
            for (int colIdx = 0; colIdx < colNumber; colIdx++)
            {
                localResult[rowIdx] += localMatrix[rowIdx * colNumber + colIdx] * vectorX[colIdx];
            }
        }

        MPI_Gather(localResult, rowsPerProcess, MPI_DOUBLE,
                resultVector, rowsPerProcess, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

        end_time = MPI_Wtime();
        elapsed_time = end_time - start_time;

        if (processRank == 0)
        {
            printf("Tempo decorrido: %f segundos\n", elapsed_time);

            printf("Cálculo executado com %d linhas, %d colunas e %d processos.\n",
                rowNumber, colNumber, numberOfProcesses);

            free(matrix);
            free(resultVector);
        }

        free(vectorX);
        free(localMatrix);
        free(localResult);

        MPI_Finalize();
        return 0;
    }