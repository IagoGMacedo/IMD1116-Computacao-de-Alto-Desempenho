import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'  # Resolve problema de threads

try:
    import numpy as np
except ImportError:
    raise ImportError("NumPy não pode ser carregado. Tente reinstalá-lo.")

def load_results(filename, shape=(381, 381, 381)):
    """Carrega resultados de forma mais eficiente"""
    try:
        return np.loadtxt(filename).reshape(shape)
    except Exception as e:
        print(f"Erro ao carregar {filename}: {str(e)}")
        raise

def calculate_mean_error(file1, file2):
    print("Carregando arquivos...")
    data1 = load_results(file1)
    data2 = load_results(file2)
    
    print("Calculando erros...")
    absolute_error = np.abs(data1 - data2)
    
    return {
        'mean_error': np.nanmean(absolute_error),
        'max_error': np.nanmax(absolute_error),
        'std_error': np.nanstd(absolute_error),
        'max_error_pos': np.unravel_index(np.nanargmax(absolute_error), absolute_error.shape)
    }

if __name__ == "__main__":
    try:
        metrics = calculate_mean_error("resultados_gpu.txt", "resultados_cpu.txt")
        print("\nMétricas de Erro:")
        print(f"Erro médio: {metrics['mean_error']:.3e}")
        print(f"Erro máximo: {metrics['max_error']:.3e} (posição: {metrics['max_error_pos']})")
        print(f"Desvio padrão: {metrics['std_error']:.3e}")
    except Exception as e:
        print(f"\nErro durante execução: {str(e)}")