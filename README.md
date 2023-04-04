# Data

The BayPass raw data files are too large to store on GitHub. BayPass was run with default parameters for the core and AUX model.
The data files are organised into:
1) data_raw - including a master files of the BayPass results (XtX and Bayes Factor) for all considered chromosomes (file SNP_XtX_BF_ALL_Data.txt.zip), with raw variants listed for results with Bayes Factor > 20 (files PC#_bfgt20_rawVariantsALL.txt). The genes were annotated using the reference genome in a .bed file (formatted from GFF using BEDOPS) and each variant was annoatated to the nearest gene using bedtools closest function.
2) gene score - R script and data files to produce gene scores -that take into account how many SNPs are annotated per gene (output files BayPassSums2_PC#.txt).
3) GOTerm - R script to run TopGO with the gene universe annotated for GO Terms (file GOTerms_desc_tab.txt), with gene pathways for the identified significant GO Terms(files SVgenetoGOmappingPC#.txt)
