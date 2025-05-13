#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define N 1000000000

static int caso2(int *vetor) {
    int soma = 0;
    for (int i = 0; i < N; i++) {
        soma += vetor[i];
    }
    return soma;
}

static int caso3(int *vetor) {
    int soma1 = 0, soma2 = 0, soma3 = 0, soma4 = 0;
    for (int i = 0; i < N; i += 4) {
        soma1 += vetor[i];
        soma2 += vetor[i+1];
        soma3 += vetor[i+2];
        soma4 += vetor[i+3];
    }
    return soma1 + soma2 + soma3 + soma4;
}

int main() {
    int *vetor = (int*)malloc(N * sizeof(int));
    clock_t inicio;
    double tempo;

    inicio = clock();
    for (int i = 0; i < N; i++) {
        vetor[i] = i * 2;
    } 
    tempo = ((double)(clock() - inicio)) / CLOCKS_PER_SEC;
    printf("Caso 1 (Inicialização): %.10f segundos\n", tempo);


    // Medição do Caso 2
    inicio = clock();
    int soma2 = caso2(vetor);
    tempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
    printf("Caso 2: %.10f segundos\n", tempo);

    // Medição do Caso 3
    inicio = clock();
    int soma3 = caso3(vetor);
    tempo = (double)(clock() - inicio) / CLOCKS_PER_SEC;
    printf("Caso 3: %.10f segundos\n", tempo);

    free(vetor);
    return 0;
}