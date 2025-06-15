import matplotlib.pyplot as plt
import numpy as np

# Dados fornecidos
sizes = [8, 64, 512, 4096, 32768, 262144, 1048576]
labels = ['8 B', '64 B', '512 B', '4 KB', '32 KB', '256 KB', '1 MB']
times = [0.000009486, 0.000009350, 0.000009703, 0.000014416, 0.000036409, 0.000150434, 0.000398953]

plt.figure(figsize=(8,5))
plt.plot(labels, times, marker='o', linestyle='-', color='b')
plt.xlabel('Tamanho da mensagem')
plt.ylabel('Tempo médio por troca (s)')
plt.title('Tempo médio de comunicação MPI vs. Tamanho da mensagem')
plt.grid(True, linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.show()