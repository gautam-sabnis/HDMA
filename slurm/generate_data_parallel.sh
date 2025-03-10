#!/bin/bash

#SBATCH -J SimData

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128GB
#SBATCH --cpus-per-task=30
#SBATCH --time=00:10:00
#SBATCH --output=logs/output_%j.out
#SBATCH --error=logs/error_%j.err
#SBATCH --array=1-300
module load apptainer

singularity exec /projects/kumar-lab/sabnig/HDMA/containers/RBayes.sif \
Rscript /projects/kumar-lab/sabnig/HDMA/mediation_DNAm/simulation_scripts/generate_data_parallel2.R  -rng 123 -P 2000 -N 2500 -rep 100 -id ${SLURM_ARRAY_TASK_ID}



