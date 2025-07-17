#!/bin/bash
#SBATCH --job-name=heat512_1n2p
#SBATCH --output=saida_%j.out
#SBATCH --error=erro_%j.err
#SBATCH --partition=intel-128
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
N=512

mpirun -np $SLURM_NTASKS ./heat2d $N