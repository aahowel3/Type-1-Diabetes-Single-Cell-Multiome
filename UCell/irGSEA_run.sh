#!/bin/bash
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBACTH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=7-00:00:00
#SBATCH --account=csd854
#SBATCH --mem=200GB

module load singularitypro
module load cpu gcc gsl

#this HAS to be run in your rstudio_server folder bc that's where all your packages are install
singularity exec --bind /cm/shared,/tscc/projects,/tscc/lustre,run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf /cm/shared/apps/containers/test/rstudio-4.3.2.sif Rscript calculate_AUCell.R
