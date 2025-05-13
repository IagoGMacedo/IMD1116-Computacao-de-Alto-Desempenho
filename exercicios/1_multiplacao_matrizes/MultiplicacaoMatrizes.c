#include <stdlib.h>
#include <time.h>
#include <stdio.h>

int **allocate_matrix(int rows, int cols);
void fill_matrix_random(int **matrix, int rows, int cols);
void print_matrix(int **matrix, int rows, int cols);
void free_matrix(int **matrix, int rows);
void multiply_with_column_access(int **A, int **B, int **C, int rowsA, int colsA, int colsB);
void multiply_with_row_access(int **A, int **B, int **C, int rowsA, int colsA, int colsB);

int main()
{
    srand(time(NULL));

    int rowsA = 1500;
    int colsA = 1500;

    int rowsB = colsA;
    int colsB = 1500;

    int **A = allocate_matrix(rowsA, colsA);
    int **B = allocate_matrix(rowsB, colsB);

    int **C_linha = allocate_matrix(rowsA, colsB);
    int **C_coluna = allocate_matrix(rowsA, colsB);

    fill_matrix_random(A, rowsA, colsA);
    fill_matrix_random(B, rowsB, colsB);

    // printf("Matriz A:\n");
    // print_matrix(A, rowsA, colsA);

    // printf("Matriz B:\n");
    // print_matrix(B, rowsB, colsB);

    multiply_with_row_access(A, B, C_linha, rowsA, colsA, colsB);

    multiply_with_column_access(A, B, C_coluna, rowsA, colsA, colsB);

    // printf("Matriz C:\n");
    // print_matrix(C, rowsA, colsB);

    // Liberar memória
    free_matrix(A, rowsA);
    free_matrix(B, rowsB);
    free_matrix(C_linha, rowsA);
    free_matrix(C_coluna, rowsA);

    return 0;
}

int **allocate_matrix(int rows, int cols)
{
    int **matrix = (int **)malloc(rows * sizeof(int *));
    if (matrix == NULL)
    {
        return NULL;
    }

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (int *)malloc(cols * sizeof(int));
        if (matrix[i] == NULL)
        {
            for (int j = 0; j < i; j++)
            {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }

    return matrix;
}

void fill_matrix_random(int **matrix, int rows, int cols)
{
    int min = 0;
    int max = 100;

    if (matrix == NULL)
        return;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            matrix[i][j] = rand() % (max - min + 1) + min;
        }
    }
}

void print_matrix(int **matrix, int rows, int cols)
{
    if (matrix == NULL)
    {
        printf("Matriz inválida!\n");
        return;
    }

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%d\t", matrix[i][j]);
        }
        printf("\n");
    }
}

void free_matrix(int **matrix, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void multiply_with_row_access(int **A, int **B, int **C, int rowsA, int colsA, int colsB)
{
    clock_t start = clock();
    for (int i = 0; i < rowsA; i++)
    {
        for (int j = 0; j < colsB; j++)
        {
            C[i][j] = 0;
            for (int k = 0; k < colsA; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    clock_t end = clock();
    double tempo_linha = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tempo método linha: %.6f segundos\n", tempo_linha);
}

void multiply_with_column_access(int **A, int **B, int **C, int rowsA, int colsA, int colsB)
{
    clock_t start = clock();
    for (int j = 0; j < colsB; j++)
    {
        for (int i = 0; i < rowsA; i++)
        {
            C[i][j] = 0;
            for (int k = 0; k < colsA; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    clock_t end = clock();
    double tempo_linha = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tempo método coluna: %.6f segundos\n", tempo_linha);
}
