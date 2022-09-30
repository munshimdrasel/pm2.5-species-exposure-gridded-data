#!/bin/bash

#SBATCH --job-name=pm.species.ec
#SBATCH --partition=normal
#SBATCH --output=/scratch/%u/%x-%N-%j.out  # Output file
#SBATCH --error=/scratch/%u/%x-%N-%j.err   # Error file
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mrasel@gmu.edu
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=32
#SBATCH --nodes=2

cd /projects/HAQ_LAB/mrasel/R/pm2.5-species-exposure-gridded-data/R


module load r-disperseR/0.1.0  

Rscript --no-restore --quiet --no-save ec.R


