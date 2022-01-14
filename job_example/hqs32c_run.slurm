#!/bin/bash

#SBATCH --job-name hqs_32c
#SBATCH --error hqs_32c_err_%j.txt
#SBATCH --output hqs_32c_out_%j.txt

#SBATCH --mail-type ALL
#SBATCH --mail-user yourmail@yourdomain

#SBATCH --partition allgroups

#SBATCH --ntasks 32
#SBATCH -c 1
#SBATCH --mem 32G
#SBATCH --time 30:00

spack load intel-parallel-studio@professional.2019.4

mpirun -np 32 ../../hqs.o random_50.txt scores_32_O0.csv -save



