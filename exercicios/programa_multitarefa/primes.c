#include <stdio.h>
#include <math.h>
#include <sys/time.h>

// gcc -DMAX=1000 -o primes primes.c -lm
// gcc -DMAX=100000000 -o primes primes.c -lm -fopenmp

// Função para verificar se um número é primo
int is_prime(int n)
{
    if (n <= 1)
        return 0;
    if (n == 2)
        return 1;
    if (n % 2 == 0)
        return 0;

    for (int i = 3; i <= sqrt(n); i += 2)
    {
        if (n % i == 0)
            return 0;
    }
    return 1;
}

int main()
{
#ifndef MAX
    printf("Erro: Defina o valor máximo com -DMAX=valor durante a compilação!\n");
    return 1;
#endif

    struct timeval start, end;
    gettimeofday(&start, NULL); // Marca o início do cálculo

    int count = 0;
    printf("Calculando primos de 2 até %d...\n", MAX);

#pragma omp parallel for
    for (int i = 2; i <= MAX; i++)
    {
        if (is_prime(i))
        {
            count++;
        }
    }

    gettimeofday(&end, NULL); // Marca o fim do cálculo

    double elapsed = (end.tv_sec - start.tv_sec) +
                     (end.tv_usec - start.tv_usec) / 1000000.0;

    printf("Total de números primos: %d\n", count);
    printf("Tempo de processamento: %.3f segundos\n", elapsed);
    return 0;
}