sequencial: gcc Main1.c -o simula -lm

paralelizado: gcc Main2.c -o simula -lm -fopenmp

exibição: python3 Exibicao.py