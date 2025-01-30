#!/bin/bash
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBACTH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --account=csd854
#SBATCH --mem=1000GB
#SBATCH --array=1-14

#cd in directory /tscc/projects/ps-gaultonlab/abhowell/scenic_data2
#conda activate pycistopic_polars
config=/tscc/projects/ps-gaultonlab/abhowell/scenic_data2/pycistopic_qc_config_rerun.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

pycistopic qc run --fragments /tscc/projects/ps-gaultonlab/abhowell/scenic_data/pipeline_out/012824_8lanes_merged.atac_fragments_new.tsv.gz --regions outs/consensus_peak_calling/consensus_regions.bed --tss qc/tss.bed --output /tscc/projects/ps-gaultonlab/abhowell/scenic_data2/qc/${sample}

