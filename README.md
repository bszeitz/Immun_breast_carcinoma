# Immune_breast_carcinoma

This is the custom R code used for the NanoString data analysis in:

Szeitz, B.; Pipek, O.; Kulka, J.; Szundi, C.; Rusz, O.; Tőkés, T.; Szász, A.M.; Kovács, K.A.; Pesti, A.; Ben Arie, T.B.; Gángó, A.; Fülöp, Z.; Drágus, E.; Vári-Kakas, S.A.; Tőkés, A.M. Investigating the Prognostic Relevance of Tumor Immune Microenvironment and Immune Gene Assembly in Breast Carcinoma Subtypes. Cancers 2022, 14, 1942. https://doi.org/10.3390/cancers14081942

**This repository consists of 2 scripts:**

| Script name  | Description  |
|---|---|
| Customized_functions.R  | This contains custom functions used in the data analysis script. |
|  NanoString_Data_Analyis_Script.R | This contains the all data analysis steps.  |

**Data analysis steps:**
1) Technical QC of the samples [1]
2) Selection of housekeeping genes for normalization [1]
3) Normalization using RUVSEq method [1, 2, 3, 4]
5) Assessment of normalization
6) Differential expression analysis with DESeq2 [4]
7) Overrepresentation analysis for significant genes
8) Volcano plots for the differential expression results
9) Unsupervised clustering of the gene expression matrix
10) Overrepresentation analysis for gene clusters
11) Export results

**References:**

[1] Bhattacharya A, Hamilton AM, Furberg H, Pietzak E, Purdue MP, Troester MA, Hoadley KA, Love MI: An approach for normalization and quality control for NanoString RNA expression data. Brief Bioinform 2021, 22(3).

[2] Bullard JH, Purdom E, Hansen KD, Dudoit S: Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments. BMC Bioinformatics 2010, 11:94.

[3] Risso D, Ngai J, Speed TP, Dudoit S: Normalization of RNA-seq data using factor analysis of control genes or samples. Nat Biotechnol 2014, 32(9):896-902.

[4] Love MI, Huber W, Anders S: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 2014, 15(12):550.

