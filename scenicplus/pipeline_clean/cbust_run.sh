#!/usr/bin/env bash

#SBATCH --job-name=cbust
#SBATCH --account=csd854
#SBATCH --qos=condo
#SBATCH --partition=condo
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time=7-00:00:00


#ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt
#mv motifs.txt cistarget_genomes/
#conda activate scenicplus
#conda install -c conda-forge python-flatbuffers
#have to modify line 301 to cluster_buster_path
#the -o argument ends in the PREFIX you want to give filenames not final folder loc
/tscc/nfs/home/abhowell/create_cisTarget_databases/create_cistarget_motif_databases.py -f /tscc/projects/ps-gaultonlab/projects/npod1/cellranger_output/Multiome/Combined/outs/Combined.with_1kb_bg_padding.fa -M /tscc/nfs/home/abhowell/aertslab_motif_colleciton/v10nr_clust_public/singletons -m /tscc/nfs/home/abhowell/cistarget_genomes/motifs.txt -o /tscc/projects/ps-gaultonlab/projects/npod1/cellranger_output/Multiome/Combined/outs/Combined --bgpadding 1000 -t 40
