#!/bin/bash
#SBATCH --job-name=cluster_sig_5.5_epoch_100
#SBATCH --output=quanah_log/%x.o%j
#SBATCH --error=quanah_log/%x.e%j
#SBATCH --partition quanah
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=5
#SBATCH --mail-user=<madison.howard@ttu.edu>
#SBATCH --mail-type=ALL

# actual code to execute
./run.sh

