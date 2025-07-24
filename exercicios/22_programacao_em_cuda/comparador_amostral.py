import random
import os
import time

def load_data(filename, expected_size):
    """Carrega os dados do arquivo como uma lista de floats."""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            data.append(float(line.strip()))
    if len(data) != expected_size:
        raise ValueError(f"Tamanho incorreto: esperado {expected_size}, obtido {len(data)}")
    return data

def compare_samples(file1, file2, shape=(381, 381, 381), sample_frac=0.01, threshold=1e-10):
    total_points = shape[0] * shape[1] * shape[2]
    sample_size = max(1000, int(total_points * sample_frac))
    
    print(f"Carregando arquivos...")
    data1 = load_data(file1, total_points)
    data2 = load_data(file2, total_points)
    
    print(f"Comparando {sample_size} pontos...")
    errors = []
    significant_diffs = []
    
    for _ in range(sample_size):
        idx = random.randint(0, total_points - 1)
        diff = abs(data1[idx] - data2[idx])
        errors.append(diff)
        if diff > threshold:
            z = idx // (shape[1] * shape[0])
            y = (idx % (shape[1] * shape[0])) // shape[0]
            x = idx % shape[0]
            significant_diffs.append({
                'position': (z, y, x),
                'error': diff
            })
    
    mean_error = sum(errors) / len(errors)
    max_error = max(errors)
    
    print("\n=== Resultados ===")
    print(f"Erro médio: {mean_error:.3e}")
    print(f"Erro máximo: {max_error:.3e}")
    print(f"Diferenças significativas (> {threshold:.1e}): {len(significant_diffs)}")
    
    if significant_diffs:
        print("\nExemplos de diferenças:")
        for diff in significant_diffs[:5]:
            print(f"Posição {diff['position']}: Erro = {diff['error']:.3e}")

if __name__ == "__main__":
    compare_samples(
        file1="resultados_gpu.txt",
        file2="resultados_cpu.txt",
        shape=(381, 381, 381),
        sample_frac=0.01
    )