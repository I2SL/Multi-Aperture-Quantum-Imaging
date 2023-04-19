#!/bin/bash
### The account associated with our lab
#SBATCH --account=ashoka
### use standard partition to get priority. Otherwise use windfall
#SBATCH --partition=standard
### Time allowed for job.
#SBATCH --time=00-96:00:00
### Number of nodes assigned to EACH JOB. Not the total number of nodes used.
#SBATCH --nodes=1
### number of tasks
#SBATCH --ntasks=1
### number of cpus per task
#SBATCH --cpus-per-task=94
### the name of the job
#SBATCH --job-name=MultiAp_CentroidSurvey
### memory per CPU
#SBATCH --mem-per-cpu=5gb
### total memory requested = memory per cpu x cpus per task. Better not to use if you have ntasks and cpus per task set
### SBATCH -mem = 470gb
### email the user once job is done
#SBATCH --mail-user=ndeshler@arizona.edu
### setup job array
#SBATCH --array 337,349,205,493,389,317,233,117,237,146,74,326,174,366,198,54,471,443,364,124,148,328,272,128,140,152,488,192,300,84
### output for a job array
#SBATCH --output=MultiAp_cfg-%a.out

#########################
### Additional flags ####
#########################
### Request GPU resources. Up to four GPUs may be requested on Puma on a single GPU node with --gres=gpu:1, 2, 3, or 4
### SBATCH --gres=gpu:1


####################################
### USEFUL ENVIRONMENT VARIABLES ###
####################################
### https://blog.ronin.cloud/slurm-job-arrays/
### $SLURM_ARRAY_JOB_ID
### $SLURM_ARRAY_TASK_ID
### $SLURM_ARRAY_TASK_MAX
### $SLURM_ARRAY_TASK_MIN


####################
### USEFUL LINKS ###
####################
### https://blog.ronin.cloud/slurm-job-arrays/



module load matlab

matlab matlab -nodisplay -nosplash -batch "Ap9Survey($SLURM_ARRAY_TASK_ID,94)" output_{$SLURM_ARRAY_TASK_ID}.txt
