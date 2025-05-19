#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main() {
    omp_set_num_threads(4);

    long long int num_pontos = 10000000;
    long long int acertos = 0;

    double inicio = omp_get_wtime();

    #pragma omp parallel
    {
        unsigned int seed = omp_get_thread_num();

        #pragma omp for reduction(+:acertos)
        for (long long int i = 0; i < num_pontos; i++) {
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;
            if (x*x + y*y <= 1.0)
                acertos++;
        }
    }

    double pi = 4.0 * acertos / num_pontos;
    double fim = omp_get_wtime();

    printf("Pi estimado: %f\n", pi);
    printf("Tempo: %f segundos\n", fim - inicio);

    return 0;
}