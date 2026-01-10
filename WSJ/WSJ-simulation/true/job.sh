#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=suntiansheng
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=24:00:00
#SBATCH --output=out/true_setting.out
#SBATCH --error=err/true_setting.err
#SBATCH --mail-user=as_448@usc.edu
#SBATCH --mail-type=END,FAIL

module purge
module load gcc/13.3.0
module load openblas/0.3.28
module load r/4.4.1

Rscript true_model_simulation.R
