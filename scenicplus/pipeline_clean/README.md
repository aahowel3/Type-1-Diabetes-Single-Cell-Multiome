R script strips out from a final annotated object the individual pieces needed going into pycistopic, notebook includes step 1 and step 2 for pycistopic and the RNA preprocessing, cbust script for step 3 to run cistarget, and the final step 4 is creating a config file and running snakemake to run SCENIC+.
For the R script - for npod1 there was an RNA object and a separate ATAC object from the same multiome data - the combined RNA+ATAC object didn't have the cell type annotations. If you have a normal multiome object ignore lines 47-68 to combine RNA + ATAC back together via shared barcodes.

Runtime 2-3 days. 
Step 1: https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html
Step 2: https://scenicplus.readthedocs.io/en/latest/human_cerebellum_scRNA_pp.html#Preprocessing-the-scRNA-seq-data
Step 3: https://scenicplus.readthedocs.io/en/latest/human_cerebellum_ctx_db.html#Creating-custom-cistarget-database
Step 4: https://scenicplus.readthedocs.io/en/latest/human_cerebellum.html#Running-SCENIC+
