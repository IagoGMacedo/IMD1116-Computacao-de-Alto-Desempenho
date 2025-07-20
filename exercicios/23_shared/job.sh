#!/bin/bash
#SBATCH --job-name=simulacao3d_cuda
#SBATCH --time=0-0:15
#SBATCH --partition=gpu-4-a100
#SBATCH --gpus-per-node=1

export TORCH_CUDA_VERSION=cu117

ulimit -s unlimited
module load compilers/nvidia/nvhpc/24.11

nsys profile -o profile_nsys ./simulacao3d_cuda