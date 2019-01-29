#!/bin/bash
#SBATCH -p general
#SBATCH -o compute.out
#SBATCH -n 50

echo -n > compute.out

module purge
module load gcc/5.4.0-alt mpi/openmpi/1.10.1-gcc r/3.5.1-openblas-gcc540

Rscript 01-16-2019_data_generation_practice.R
