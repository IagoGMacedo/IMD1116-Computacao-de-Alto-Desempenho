#include <stdio.h>
#include <math.h>
#include <time.h>

#define ITERATIONS 10000000000L

int main() {
    clock_t start = clock();

    double sum = 0.0;
    double sign = 1.0;

    for (long k = 0; k < ITERATIONS; k++) {
        sum += sign / (2 * k + 1);
        sign *= -1;
    }

    double pi_approx = sum * 4;

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    double error = fabs(pi_approx - M_PI);

    printf("Iterations: %ld\n", ITERATIONS);
    printf("Approximation: %.15f\n", pi_approx);
    printf("Real value:    %.15f\n", M_PI);
    printf("Absolute error: %.15f\n", error);
    printf("Execution time: %.6f seconds\n", time_spent);

    return 0;
}