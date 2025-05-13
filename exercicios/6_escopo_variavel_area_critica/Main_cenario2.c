#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

int main()
{
    long long int total_pontos = 1000000;
    long long int dentro_circulo = 0;

    srand(time(NULL));

    #pragma omp parallel shared(total_pontos, dentro_circulo)
    {
        long long int local_dentro = 0;

        unsigned int seed = (unsigned int)time(NULL) ^ omp_get_thread_num();

        #pragma omp for private(seed) firstprivate(total_pontos) lastprivate(local_dentro)
        for (long long int i = 0; i < total_pontos; i++)
        {
            double x = (double)rand_r(&seed) / RAND_MAX;
            double y = (double)rand_r(&seed) / RAND_MAX;

            if ((x * x + y * y) <= 1.0)
            {
                local_dentro++;
            }
        }

        #pragma omp critical
        {
            dentro_circulo += local_dentro;
        }
    }

    double pi_estimado = 4.0 * dentro_circulo / total_pontos;
    printf("Estimativa de pi: %.10f\n", pi_estimado);

    return 0;
}