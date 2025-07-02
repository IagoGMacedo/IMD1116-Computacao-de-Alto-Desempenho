#include <stdio.h>
#include <sys/time.h>

#define TOL 1e-6

int main()
{
    int N = 100000000;

    float a[N], b[N], c[N], res[N];
    int err = 0;
    struct timeval start, end;
    double elapsed;

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        a[i] = (float)i;
        b[i] = 2.0 * (float)i;
        c[i] = 0.0;
        res[i] = i + 2 * i;
    }

    gettimeofday(&start, NULL);

#pragma omp parallel for
    for (int i = 0; i < N; i++)
    {
        c[i] = a[i] + b[i];
    }

    gettimeofday(&end, NULL);

    elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_usec - start.tv_usec) / 1000000.0);

    printf("Tempo para soma (versão com cpu): %f segundos\n", elapsed);

#pragma omp parallel for reduction(+ : err)
    for (int i = 0; i < N; i++)
    {
        float val = c[i] - res[i];
        val = val * val;
        if (val > TOL)
            err++;
    }
    printf("vectors added with %d errors\n", err);
    return 0;
}