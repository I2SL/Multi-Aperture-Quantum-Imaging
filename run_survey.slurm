#!/bin/bash
### The account associated with our lab
#SBATCH --account=ashoka
### use standard partition to get priority. Otherwise use windfall
#SBATCH --partition=standard
### Time allowed for job.
#SBATCH --time=72:00:00
### Number of nodes assigned to EACH JOB. Not the total number of nodes used.
#SBATCH --nodes=1
### number of tasks
#SBATCH --ntasks=1
###
#SBATCH --cpus-per-task=94
###
#SBATCH --job-name=MultiApSurvey
### memory per CPU
#SBATCH --mem-per-cpu=5gb
###
### SBATCH -mem = 512gb
###
###
#SBATCH --mail-user=ndeshler@arizona.edu
#SBATCH --array 1-

#################################
### USEFUL ENVIRONMENT VARIABLES
#################################
### $SLURM_ARRAY_JOB_ID
### $SLURM_ARRAY_TASK_ID
### $SLURM_ARRAY_TASK_MAX
### $SLURM_ARRAY_TASK_MIN




module load matlab

matlab -nodisplay -nosplash < MultiAperture_Survey($SLURM_ARRAY_TASK_ID) > survey_output.txt
