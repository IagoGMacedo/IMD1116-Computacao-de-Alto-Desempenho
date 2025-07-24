#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 100000000
#define TOL 0.0000001

int main() {
    float *a, *b, *c, *res;
    int err = 0;

    // Alocação dinâmica na heap
    a = (float *) malloc(N * sizeof(float));
    b = (float *) malloc(N * sizeof(float));
    c = (float *) malloc(N * sizeof(float));
    res = (float *) malloc(N * sizeof(float));

    if (a == NULL || b == NULL || c == NULL || res == NULL) {
        printf("Erro ao alocar memória.\n");
        return 1;
    }

    double init_time, compute_time, test_time;

    init_time = -omp_get_wtime();

    // Inicialização dos vetores
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        a[i] = (float)i;
        b[i] = 2.0f * (float)i;
        c[i] = 0.0f;
        res[i] = i + 2 * i;
    }

    init_time += omp_get_wtime();
    compute_time = -omp_get_wtime();

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        c[i] = a[i] + b[i];
    }

    compute_time += omp_get_wtime();
    test_time = -omp_get_wtime();

    // Verificação de erros
    #pragma omp parallel for reduction(+:err)
    for (int i = 0; i < N; i++) {
        float val = c[i] - res[i];
        if (val * val > TOL)
            err++;
    }

    test_time += omp_get_wtime();

    printf("Vectors added with %d errors\n", err);
    printf("Init time:    %.3fs\n", init_time);
    printf("Compute time: %.3fs\n", compute_time);
    printf("Test time:    %.3fs\n", test_time);
    printf("Total time:   %.3fs\n", init_time + compute_time + test_time);

    // Liberação da memória
    free(a);
    free(b);
    free(c);
    free(res);

    return 0;
}
