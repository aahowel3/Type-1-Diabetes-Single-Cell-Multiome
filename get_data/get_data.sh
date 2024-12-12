#!/bin/bash
#SBATCH --partition=condo
#SBATCH --qos=condo
#SBACTH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00:00
#SBATCH --account=csd854
#SBATCH --array=1-2

#mamba init
#mamba activate get_data
#the nodes and cpus are for each job rather than the whole array
#config=/tscc/nfs/home/abhowell/get_data/srrs_rnaseqonly_crossref.config

config=/tscc/nfs/home/abhowell/get_data/srrs_ATAConly_crossref_sizes.config
#try temp config where gse that has 100GB fastqs arent in list
#config=/tscc/nfs/home/abhowell/get_data/temp.txt
SRX=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$3==ArrayTaskID {print $1}' $config)
SRR=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$3==ArrayTaskID {print $2}' $config)

mkdir -p /tscc/projects/ps-gaultonlab/abhowell/scRNA-seq/${SRX}
cd /tscc/projects/ps-gaultonlab/abhowell/scRNA-seq/${SRX}
prefetch $SRR --max-size 420000000000
parallel-fastq-dump --sra-id $SRR --threads 8 --outdir /tscc/projects/ps-gaultonlab/abhowell/scRNA-seq/${SRX} --split-files --gzip --tmpdir /tscc/projects/ps-gaultonlab/abhowell/scRNA-seq/${SRX}

