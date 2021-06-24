#!/bin/bash
#SBATCH --job-name=RunRun_matador
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40

# actual code to execute
./run.sh
