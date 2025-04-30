#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

int main()
{
    long long int total_pontos = 1000000;
    long long int dentro_circulo = 0;

    srand(time(NULL));

    #pragma omp parallel for
    for (long long int i = 0; i < total_pontos; i++)
    {
        double x = (double)rand() / RAND_MAX;
        double y = (double)rand() / RAND_MAX;

        if ((x * x + y * y) <= 1.0)
        {
            dentro_circulo++;
        }
    }

    double pi_estimado = 4.0 * dentro_circulo / total_pontos;
    printf("Estimativa de pi: %.10f\n", pi_estimado);

    return 0;
}