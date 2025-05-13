#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// gcc -fopenmp -O2 Main2.c -o pi_critical

int main() {
    omp_set_num_threads(16);

    long long int num_pontos = 10000000;
    int num_threads;
    long long int *acertos_thread;

    double inicio = omp_get_wtime();

    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        unsigned int seed = id;

        #pragma omp single
        {
            num_threads = omp_get_num_threads();
            acertos_thread = calloc(num_threads, sizeof(long long int));
        }

        long long int acertos_local = 0;

        #pragma omp for
        for (long long int i = 0; i < num_pontos; i++) {

            // com rand
            // double x = (double)rand() / RAND_MAX;
            // double y = (double)rand() / RAND_MAX;

            // com rand_r
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;

            if (x*x + y*y <= 1.0)
                acertos_local++;
        }

        acertos_thread[id] = acertos_local;
    }

    long long int acertos = 0;
    for (int i = 0; i < num_threads; i++)
        acertos += acertos_thread[i];

    double pi = 4.0 * acertos / num_pontos;
    double fim = omp_get_wtime();

    printf("Pi estimado: %f\n", pi);
    printf("Tempo: %f segundos\n", fim - inicio);

    free(acertos_thread);
    return 0;
}