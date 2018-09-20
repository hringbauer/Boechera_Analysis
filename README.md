# Boechera_Analysis
Analysis of Boechera Stricta Data

The raw and processed data is found in various folders in /Data:

***IBD_raw***:
Contains the raw .ibd files (i.e. Output from Beagle)

***LinkageMap***:
Contains the linkage map of Boechera Stricta

***IBD_cM***:
Contains .ibd files transformed to cM based on linkage map

***IBD_merged***:
Contains .ibd files that have been merged to bridge spurious gaps, as well as overlaps.


The **Snakefile** in main runs the Python scripts for these two jobs (transforming to cM, as well as merging gaps) when needed.



*For Analysis of the process Data*:
There are 2 Jupyter Notebooks in Notebooks:

1) **first_analysis.ipynb**
Load the data (blocks sharing and individuals coordinates), bring it in a form useable for Python and do some preliminary analysis

2) **snp_freqs.ipynb**
Load the unphased, row SNP-data. Do preliminary analysis of SNPs: Site Frequency Spectrum, Heterozygosity to estimate fundamental patterns such as F_ST. Use scikit-allele to estimate F_IS using Weir-Cockerhams estimator.

In **Image_Analysis**:
Some figures of the output



All rights reserved

