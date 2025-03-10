#!/bin/bash

#SBATCH -J BOOM

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=512GB
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00
#SBATCH --output=logs/output_%j.out
#SBATCH --error=logs/error_%j.err
#SBATCH --array=1-300

module load apptainer

singularity exec /projects/kumar-lab/sabnig/HDMA/containers/RBayes.sif \
Rscript /projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_scripts/boom_parallel.R  -rep 100 -id ${SLURM_ARRAY_TASK_ID}



