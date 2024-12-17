#installed cwtool into base mamba (shouldnt have done that probs)
#docker is already available and ready to run
#docker WILL NOT run on TSCC
git clone https://github.com/hubmapconsortium/salmon-rnaseq.git
#requires a checked out version to run so docker images are up to date
tmux attach -t hubmap 
cd salmon_rnaseq 
git checkout v2.0.6
cwltool pipeline.cwl --assay 10x_v3 --fastq_dir /nfs/lab/projects/mega_pancreas/data/fastq/HPAP/RNA/HPAP-101 --threads 1
