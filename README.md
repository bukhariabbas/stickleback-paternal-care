# stickleback-paternal-care

This contains R scripts which were used to detect differentially expressed genes using edgeR. The files of Rscripts either contain SO or SC. SO stands for social opportunity (Stickleback paternal care), whereas SC stands for social challenge (stickleback territorial challenge). 

For multiple testing correction, we implemented empirical FDR. R scripts with "permNulldist" tags were used to generate the null distribution of each contrast, whereas R script with "perm.Pval" tags were used to compute an FDR.
