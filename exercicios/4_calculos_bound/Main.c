#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>  // Adicionada a biblioteca OpenMP

//rodando:
// gcc -o main Main.c -lm -fopenmp
// ./main 100

// Função para alocar matriz
int **allocate_matrix(int rows, int cols) {
    int **matrix = (int **)malloc(rows * sizeof(int *));
    for (int i = 0; i < rows; i++) {
        matrix[i] = (int *)malloc(cols * sizeof(int));
    }
    return matrix;
}

// Função para preencher matriz com valores aleatórios
void fill_matrix_random(int **matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = rand() % 100;
        }
    }
}

// Função para liberar matriz
void free_matrix(int **matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

// Função para somar as matrizes A, B e C em D
void sum_matrix(int **A, int **B, int **C, int **D, int rows, int cols) {
    #pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            D[i][j] = A[i][j] + B[i][j] + C[i][j];
        }
    }
}

double time_diff(struct timeval start, struct timeval end) {
    return (end.tv_sec - start.tv_sec) + 1e-6 * (end.tv_usec - start.tv_usec);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Uso: %s <numero_de_threads>\n", argv[0]);
        return 1;
    }

    int operationNumbers = 100000000;
    
    int num_threads = atoi(argv[1]);
    omp_set_num_threads(num_threads);  

    struct timeval start, end;
    double cpu_time, memory_time;

    srand(time(NULL));

    // Medição do tempo CPU Bound
    gettimeofday(&start, NULL);
    
    double sumt = 0.0;
    #pragma omp parallel for 
    for (int i = 0; i < operationNumbers; i++) {
        sumt += tan(sin(i) * cos(i)) + sqrt(i);
    }
    
    gettimeofday(&end, NULL);
    cpu_time = time_diff(start, end);

    int rows = 1;
    int cols = operationNumbers;

    int **A = allocate_matrix(rows, cols);
    int **B = allocate_matrix(rows, cols);
    int **C = allocate_matrix(rows, cols);
    int **D = allocate_matrix(rows, cols);

    fill_matrix_random(A, rows, cols);
    fill_matrix_random(B, rows, cols);
    fill_matrix_random(C, rows, cols);

    // Medição do tempo Memory Bound
    gettimeofday(&start, NULL);
    
    sum_matrix(A, B, C, D, rows, cols);
    
    gettimeofday(&end, NULL);
    memory_time = time_diff(start, end);

    free_matrix(A, rows);
    free_matrix(B, rows);
    free_matrix(C, rows);
    free_matrix(D, rows);

    printf("Threads usadas: %d\n", num_threads);
    printf("Tempo CPU Bound: %.4f segundos\n", cpu_time);
    printf("Tempo Memory Bound: %.4f segundos\n", memory_time);

    return 0;
}