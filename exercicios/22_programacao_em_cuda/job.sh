#!/bin/bash
#SBATCH --job-name=simulacao3d_cuda
#SBATCH --time=0-0:15
#SBATCH --partition=gpu-8-v100
#SBATCH --gpus-per-node=1

export TORCH_CUDA_VERSION=cu117

module load compilers/nvidia/nvhpc/24.11

./simulacao3d_cuda