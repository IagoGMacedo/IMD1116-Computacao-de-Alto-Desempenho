#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// gcc -fopenmp -O2 Main1.c -o pi_critical

int main() {
    omp_set_num_threads(4);

    long long int num_pontos = 10000000;
    long long int acertos = 0;

    double inicio = omp_get_wtime();

    #pragma omp parallel
    {
        unsigned int seed = omp_get_thread_num();
        long long int acertos_thread = 0;

        #pragma omp for
        for (long long int i = 0; i < num_pontos; i++) {

            // com rand
            // double x = (double)rand() / RAND_MAX;
            // double y = (double)rand() / RAND_MAX;


            // com rand_r
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;
            if (x*x + y*y <= 1.0)
                acertos_thread++;
        }

        #pragma omp critical
        acertos += acertos_thread;
    }

    double pi = 4.0 * acertos / num_pontos;
    double fim = omp_get_wtime();

    printf("Pi estimado: %f\n", pi);
    printf("Tempo: %f segundos\n", fim - inicio);

    return 0;
}