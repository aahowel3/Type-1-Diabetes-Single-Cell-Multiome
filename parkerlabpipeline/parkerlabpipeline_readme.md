#in a tmux window first 
tmux attach -t parkerlab
#everytime you log back in reset enviornment 
srun -N 1 -n 16 -t 1-00:00:00 -p condo -q condo -A csd854 --pty /bin/bash
conda activate parkerlab
module load singularitypro/3.11
#you have to run it on /projects - tmp work folder created in snRNAseq-Nextflow is too big
cd /tscc/projects/ps-gaultonlab/abhowell/snRNAseq-NextFlow

#install nextflow via mamba
conda create -n parkerlab
conda activate parkerlab
conda install bioconda::nextflow

#if you need to reclone the parkerlab repo re-copy ruths edited config files
git clone https://github.com/ParkerLab/snRNAseq-NextFlow.git
cd snRNAseq-NextFlow
cp /tscc/projects/ps-gaultonlab/relgamal/snRNAseq-NextFlow-master/library-config_ruth.json . 
cp /tscc/projects/ps-gaultonlab/relgamal/snRNAseq-NextFlow-master/3M-february-2018.txt . 
#changes are made to the slurm accout name and the hardcoded paths of the reference genome/star indicies
#location for singularity-cache is specified with all downloaded .img files
cp /tscc/projects/ps-gaultonlab/relgamal/snRNAseq-NextFlow-master/nextflow.config . 

#actual command
nextflow run -resume -params-file library-config_ruth.json --barcode-whitelist /tscc/projects/ps-gaultonlab/abhowell/snRNAseq-NextFlow/3M-february-2018.txt \
--chemistry V3 --results /tscc/projects/ps-gaultonlab/abhowell/ main.nf
