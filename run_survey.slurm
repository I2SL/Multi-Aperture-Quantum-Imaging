#!/bin/bash
#SBATCH --account=ashoka
#SBATCH --partition=standard
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --job-name=MultiApSurvey
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-user=ndeshler@arizona.edu


module load matlab

matlab -nodisplay -nosplash < MultiAperture_Survey.m > output.txt
