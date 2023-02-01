#!/bin/bash
#SBATCH --account=ashoka
#SBATCH --standard
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --job-name=MultiApSurvey
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-user=ndeshler@arizona.edu

module load matlab

export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXX)

matlab -nodisplay -nosplash run.m
