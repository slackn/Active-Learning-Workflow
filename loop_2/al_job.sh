#!/bin/bash
#SBATCH --time=48:00:00      # adjust as needed
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=10          
#SBATCH --cpus-per-task=1    # only one thread
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@sabanciuniv.edu
#SBATCH --job-name=al_Na8

module load chem/turbomole/7.9

eval "$(conda shell.bash hook)"
conda activate gpaw-env

# test Turbomole availability
which define

python active_learning_loop.py -c config.yaml
