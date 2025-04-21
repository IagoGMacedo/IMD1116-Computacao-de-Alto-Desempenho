#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Função para somar as matrizes A, B e C em D
void sum_matrix(int **A, int **B, int **C, int **D, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            D[i][j] = A[i][j] + B[i][j] + C[i][j];
        }
    }
}

int main()
{
    // CPU BOUND
    double sumt = 0.0;
    for (int i = 0; i < 1000; i++)
    {
        sumt += tan(sin(i) * cos(i)) + sqrt(i);
    }

    // MEMORY BOUND

    //alocação das matrizes
    int rows = 1500;
    int cols = 1500;

    int **A = allocate_matrix(rows, cols);
    int **B = allocate_matrix(rows, cols);
    int **C = allocate_matrix(rows, cols);
    int **D = allocate_matrix(rows, cols);


    fill_matrix_random(A, rows, cols);
    fill_matrix_random(B, rows, cols);
    fill_matrix_random(C, rows, cols);

    //soma das matrizes
    sum_matrix(A,B,C,D, rows, cols);

    //liberando matrizes
    free_matrix(A, rows);
    free_matrix(B, rows);
    free_matrix(C, rows);
    free_matrix(D, rows);

    return 0;
}
