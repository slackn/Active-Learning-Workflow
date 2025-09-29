#!/bin/bash
#SBATCH --time=10:00:00      # adjust as needed
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1          # only one task
#SBATCH --cpus-per-task=1    # only one thread
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sila.horozoglu@sabanciuniv.edu
#SBATCH --job-name=relabel

module load chem/turbomole/7.9

eval "$(conda shell.bash hook)"
conda activate gpaw-env

python relabel.py
