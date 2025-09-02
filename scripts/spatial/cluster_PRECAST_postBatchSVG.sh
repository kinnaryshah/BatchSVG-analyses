#!/bin/bash
#SBATCH --job-name=precast
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=END
#SBATCH --mail-user=kinnaryshahh@gmail.com 

echo "**** Job starts ****"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5.x

## List current modules for reproducibility
module list

## Edit with your job command
Rscript cluster_PRECAST_postBatchSVG.R

echo "**** Job ends ****"
date