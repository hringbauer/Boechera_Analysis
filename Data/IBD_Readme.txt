*** For file CB10k.ibd30.ibd.txt:

Here attached IBD of 108 CB genotypes. For each line in the file:
Column1: First sample’s name
Column2: First sample's haplotype index; 1 or 2
Column3: Second sample’s name
Column4: Second sample’s haplotype index; 1 or 2
Column5: Scaffold
Column6: Start position of the IBD
Column7: End position of the IBD
Column8: LOD score, larger values indicate greater evidence for IBD
Column8: length of the IBD

Methods: 
We imputed missing data and phased genotypes using Beagle v4.1 with default settings. After that we performed identify-by-descent (IBD) analysis using Beagle with ibdtrim = 30. For other parameters, the default settings were used. 

There are 128049 shared haplotype pairs, with median length ~ 179kb. 


*** For file CB108.impute.geno.gz:
 
In the genotype file, each line is for one SNP, the first 4 columns are: Scaffold, Position, Reference allele, Alternative allele; after that are genotypes for the 108 CB individuals.
Genotype (0, 1 and 2) coded as the number of alternative allele:
0: homozygote of reference allele
1: Heterozygote
2:homozygote of alternative allele
 
There are ~237K SNPs segregating in this population. 

*** For file CB108.coordinate.txt:

Three columns have information on ID_number, Latitude, and Longitude. 








 
