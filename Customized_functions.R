##############################
#
# Utility functions
#
# Collected by Beáta Szeitz
#
##############################


##############################
# Extract flag for binding density
# 
# Arguments for this function:
# - rcc.file: the RCC file
# - lower: the lower limit for binding density
# - upper: the upper limit for binding density
#
# This function returns a character indicating
# whether there is a flag for binding density
#
# Packages required: none
# 
binding_density_flag <- function(rcc.file, lower, upper){
  bin.dens <- as.numeric(rcc.file$Lane_Attributes[["BindingDensity"]])
  if (bin.dens < upper & bin.dens > lower){
    res <- "No flag"
  } else {
    res <- "Flag"
  }
  return(res)
}


##############################
# Extract flag for image quality
# 
# Arguments for this function:
# - rcc.file: the RCC file
#
# This function returns a character indicating
# whether there is a flag for image quality
#
# Packages required: none
# 
image_quality_flag <- function(rcc.file){
  FOV.ratio <- as.numeric(rcc.file$Lane_Attributes[["FovCounted"]]) / as.numeric(rcc.file$Lane_Attributes[["FovCount"]])
  if (FOV.ratio > 0.75){
    res <- "No flag"
  } else {
    res <- "Flag"
  }
  return(res)
}

##############################
# Extract flags for limit of detection
# 
# Arguments for this function:
# - rcc.file: the RCC file
#
# This function returns a character indicating
# whether there is a flag for limit of detection
#
# Packages required: none
# 
LoD_flag <- function(rcc.file){
  gene.counts <- rcc.file$Code_Summary
  POS_E <- as.numeric(gene.counts[gene.counts$Name =="POS_E","Count"])
  NEGs <- as.numeric(gene.counts[grep("NEG", gene.counts$Name),"Count"])
  NEGs.mean.SD <- mean(NEGs) +sd(NEGs)
  if (POS_E > NEGs.mean.SD){
    res <- "No flag"
  } else {
    res <- "Flag"
  }
  return(res)
}

##############################
# Extract flags for positive control linearity
# 
# Arguments for this function:
# - rcc.file: the RCC file
#
# This function returns a character indicating
# whether there is a flag for positive control linearity
#
# Packages required: none
# 
pos_ctrl_lin_flag <- function(rcc.file){
  gene.counts <- rcc.file$Code_Summary
  pos.ctrls <- as.numeric(gene.counts[grep("POS_", gene.counts$Name),"Count"])
  names(pos.ctrls) <- gene.counts[grep("POS_", gene.counts$Name),"Name"]
  pos.ctrls <- pos.ctrls[order(names(pos.ctrls))]
  reference <- c(128, 128/4, 128/16, 128/64, 128/256, 128/1024)
  r.sq <- summary(lm(pos.ctrls ~ reference))$r.squared
  if (r.sq > 0.95){
    res <- "No flag"
  } else {
    res <- "Flag"
  }
  return(res)
}



##############################
# Perform RUVSeq normalization
# according to the guidelines provided in 
# https://doi.org/10.1093/bib/bbaa163 and 
# https://github.com/bhattacharya-a-bt/CBCS_normalization/blob/master/CBCS_normalization_tutorial.pdf
# 
# Arguments for this function:
# - Expr: the expression table (genes in rows and samples in columns)
# - Clin: the clinical table
# - Genes: the metadata of the genes
# - k: how many dimensions of unwanted variation should be estimated
# - exclude: which housekeeping genes need to be excluded
#
# This function returns a SeqExpressionSet and the 
# variance-stabilized 
#
# Packages required: EDASeq, RUVSeq, DESeq2, limma
# 
RUVSeq_normalization <- function(Expr, Clin, Genes, k, exclude = NULL){
  
  # Ensure that gene names are in the same order
  Genes <- Genes[rownames(Expr),] 
  
  # Create newSeqExpressionSet (EDASeq package)
  set <- newSeqExpressionSet(as.matrix(round(Expr)),
                             phenoData=Clin,
                             featureData=Genes)
  
  # Specify the housekeeping genes used for normalization
  housek.used <- Genes[Genes$CodeClass =="Housekeeping","Name"]
  housek.used <- housek.used[! housek.used %in% exclude]
  
  # Upper quartile normalization followed by RUVg
  set <- betweenLaneNormalization(set, which="upper")
  set <- RUVg(set, housek.used, k=k)
  deseq.dataset <- DESeqDataSetFromMatrix(counts(set), colData=pData(set), design=~1)
  rowData(deseq.dataset) <- Genes
  
  # Estimate size factors
  deseq.dataset <- estimateSizeFactors(deseq.dataset)
  deseq.dataset <- estimateDispersionsGeneEst(deseq.dataset)
  deseq.counts <- counts(deseq.dataset, normalized=TRUE)
  disp.est <- pmax((rowVars(deseq.counts) - rowMeans(deseq.counts)),0) / rowMeans(deseq.counts)^2 
  mcols(deseq.dataset)$dispGeneEst <- disp.est
  deseq.dataset <- estimateDispersionsFit(deseq.dataset, fitType="mean")
  
  # Variance stabilizing transformation
  var.stab <- varianceStabilizingTransformation(deseq.dataset, blind=FALSE)
  expr.mat <- assay(var.stab)
  
  # Remove the variation estimated by RUVg
  est.var <- as.matrix(colData(deseq.dataset)[,grep("W_",colnames(colData(deseq.dataset)))]) # this can potentially introduce problems if any annotation has "W_" in it
  expr.mat <- removeBatchEffect(expr.mat, covariates=est.var)
  assay(var.stab) <- expr.mat
  return(list(set = set,
              vsd = var.stab))
}


##############################
# Extract DESeq2 results
# 
# Arguments for this function:
# - resultfile: a data frame containing at least 1 column containing the gene names
# - DESeqDataSet.in: created by DESeqDataSetFromMatrix function
# - design.in: the design
# - alpha.in: the significance level
# - contrast.in: the contrast to extract
# - colname.out: character string to append column names with
#
# This function returns diff. expr. results from DESEq2 in a formatted way.
#
# Packages required: DESeq2
#
extract_results_DEseq2 <- function(resultfile, DESeqDataSet.in, design.in, alpha.in, 
                                   contrast.in=NULL, colname.out){
  design(DESeqDataSet.in) <- design.in
  res <- DESeq(DESeqDataSet.in)
  res <- results(res, contrast=contrast.in, alpha = alpha.in)
  res_suppl <- cbind(res@rownames,as.data.frame(res@listData))
  colnames(res_suppl)[2:ncol(res_suppl)] <- paste(colname.out, colnames(res_suppl)[2:ncol(res_suppl)], sep="_")
  colnames(res_suppl)[1] <- "Gene"
  resultfile_out <- merge(resultfile, res_suppl, by="Gene")
  return(resultfile_out)
}


##############################
# Contingency table generation
# used for perform_fisher function within ORA_test and ORA_test_onlyBOTH functions
# 
# Arguments for this function:
# - among.sign: number of genes within the category AND are significant
# - sign.genes: number of significant genes
# - no.of.genes: number of genes within the category
#
# This function returns a contingency table used for Fisher test.
#
# Packages required: none
# 
crosstab <- function(among.sign, sign.genes,no.of.genes){
  d <- data.frame(gene.in.interest=c(among.sign, length(sign.genes)-among.sign), 
                  gene.not.interest=c(no.of.genes-among.sign,
                                      730-length(sign.genes)-no.of.genes-among.sign) 
  )
  row.names(d) <- c("In_category", "not_in_category")
  return(d)
}


##############################
# Extract Fisher's exact test results
# used for ORA_test and ORA_test_onlyBOTH functions
# 
# Arguments for this function:
# - d: the contingency table
#
# This function returns a vector of Fisher test results:
# - confidence interval, odds ratio, p-value
#
# Packages required: none
# 
perform_fisher <- function(d){
  fis <- fisher.test(d, alternative = "greater",conf.int = TRUE, conf.level = 0.95)
  CI <- paste0(round(fis$conf.int[1],4)," - ",round(fis$conf.int[2],4))
  OR <- fis[["estimate"]][["odds ratio"]]
  p <- fis$p.value
  return(c(CI, OR, p))
}


##############################
# Perform overrepresentation analysis
# for the specified genes (separately for UP/DOWN and also combined)
# 
# Arguments for this function:
# - sign.genes.UP: character vector of upregulated genes
# - sign.genes.DOWN: character vector of downregulated genes
# - Name: character string to append column names with
# - IR.categories: a matrix indicating which categories the genes are part of
#
# This function returns a data frame containing the following columns:
# - Number of genes in the category
# - Odds ratio (with 95% confidence interval)
# - Fisher test p-value
# - Fisher test adjusted p-value (method: Bonferroni)
# - The list of the differentially expressed genes present in the category
#
# Packages required: none
# 
ORA_test <- function(sign.genes.UP, sign.genes.DOWN, Name, IR.categories){
  
  sign.genes.BOTH <- c(sign.genes.UP, sign.genes.DOWN)
  panel_geneList_sub_BOTH <- IR.categories[row.names(IR.categories) %in% sign.genes.BOTH,]
  panel_geneList_sub_UP <- IR.categories[row.names(IR.categories) %in% sign.genes.UP,]
  panel_geneList_sub_DOWN <- IR.categories[row.names(IR.categories) %in% sign.genes.DOWN,]
  
  summarize_for_one <- function(sumtable){
    sumtable_sub <- sumtable
    catg_matrix <- matrix(nrow=ncol(sumtable_sub), ncol=2)
    row.names(catg_matrix) <- colnames(sumtable_sub)
    colnames(catg_matrix) <- c("Category","No.of.Genes.in.Category")
    catg_matrix[,1]<- colnames(sumtable_sub)
    
    catg_matrix[,2] <- apply(sumtable_sub, 2, function(x){length(x[x=="+"])})
    
    columns <- c("No.DEGs.in.category","OR (CI(95%))", "p.value", "p.adj", "DEGs.in.category")
    resultmat <- matrix(ncol=length(columns)*3, nrow=nrow(catg_matrix))
    colnames(resultmat)<- c(paste(columns,"BOTH", sep="_"),paste(columns,"UP", sep="_"),paste(columns,"DOWN", sep="_"))
    row.names(resultmat) <- row.names(catg_matrix)
    sumtable_sub <- as.data.frame(cbind(catg_matrix, resultmat))
    
    for (i in 1:nrow(sumtable_sub)){
      category <- sumtable_sub[i,1]
      no.of.genes <- as.numeric(sumtable_sub$No.of.Genes.in.Category[i])
      among.sign_BOTH <- nrow(panel_geneList_sub_BOTH[panel_geneList_sub_BOTH[,category]=="+",])
      sumtable_sub$DEGs.in.category_BOTH[i] <- paste(row.names(panel_geneList_sub_BOTH[panel_geneList_sub_BOTH[,category]=="+",]), collapse = "/")
      
      among.sign_UP <- nrow(panel_geneList_sub_UP[panel_geneList_sub_UP[,category]=="+",])
      sumtable_sub$DEGs.in.category_UP[i] <- paste(row.names(panel_geneList_sub_UP[panel_geneList_sub_UP[,category]=="+",]), collapse = "/")
      
      among.sign_DOWN <- nrow(panel_geneList_sub_DOWN[panel_geneList_sub_DOWN[,category]=="+",])
      sumtable_sub$DEGs.in.category_DOWN[i] <- paste(row.names(panel_geneList_sub_DOWN[panel_geneList_sub_DOWN[,category]=="+",]), collapse = "/")
      
      sumtable_sub$No.DEGs.in.category_BOTH[i] <- paste0(among.sign_BOTH, " (out of ", length(sign.genes.BOTH),")")
      sumtable_sub$No.DEGs.in.category_UP[i] <- paste0(among.sign_UP, " (out of ", length(sign.genes.BOTH),")")
      sumtable_sub$No.DEGs.in.category_DOWN[i] <- paste0(among.sign_DOWN, " (out of ", length(sign.genes.BOTH),")")
      
      if (among.sign_BOTH==0){
        sumtable_sub[i,c("OR (CI(95%))_BOTH", "p.value_BOTH")] <- NA
      } else {
        fish.res <- perform_fisher(crosstab(among.sign_BOTH, sign.genes.BOTH,no.of.genes))
        sumtable_sub[i,c("OR (CI(95%))_BOTH", "p.value_BOTH")] <- c(paste0(round(as.numeric(fish.res[2]),4)," (",fish.res[1],")"), fish.res[3])
      }
      
      if (among.sign_UP==0){
        sumtable_sub[i,c("OR (CI(95%))_UP", "p.value_UP")] <- NA
      } else {
        fish.res <- perform_fisher(crosstab(among.sign_UP, sign.genes.UP,no.of.genes))
        sumtable_sub[i,c("OR (CI(95%))_UP", "p.value_UP")] <- c(paste0(round(as.numeric(fish.res[2]),4)," (",fish.res[1],")"), fish.res[3])
      }
      
      if (among.sign_DOWN==0){
        sumtable_sub[i,c("OR (CI(95%))_DOWN", "p.value_DOWN")] <- NA
      } else {
        fish.res <- perform_fisher(crosstab(among.sign_DOWN, sign.genes.DOWN,no.of.genes))
        sumtable_sub[i,c("OR (CI(95%))_DOWN", "p.value_DOWN")] <- c(paste0(round(as.numeric(fish.res[2]),4)," (",fish.res[1],")"), fish.res[3])
      }
    }
    
    sumtable_sub$p.adj_BOTH <- p.adjust(sumtable_sub$p.value_BOTH, method="bonferroni")
    sumtable_sub$p.adj_UP <- p.adjust(sumtable_sub$p.value_UP, method="bonferroni")
    sumtable_sub$p.adj_DOWN <- p.adjust(sumtable_sub$p.value_DOWN, method="bonferroni")
    return(sumtable_sub)
  }
  category_summary_sub <- summarize_for_one(IR.categories)
  
  colnames(category_summary_sub)[3:ncol(category_summary_sub)] <- paste(Name, colnames(category_summary_sub)[3:ncol(category_summary_sub)], sep="_")
  
  return(category_summary_sub)
}


##############################
# Perform overrepresentation analysis
# for the specified genes
# 
# Arguments for this function:
# - gene.vector: character vector of the genes
# - Name: character string to append column names with
# - IR.categories: a matrix indicating which categories the genes are part of
#
# This function returns a data frame containing the following columns:
# - Number of genes in the category
# - Odds ratio (with 95% confidence interval)
# - Fisher test p-value
# - Fisher test adjusted p-value (method: Bonferroni)
# - The list of the differentially expressed genes present in the category
#
# Packages required: none
# 
ORA_test_onlyBOTH <- function(gene.vector, Name, IR.categories){
  
  panel_geneList_sub <- IR.categories[row.names(IR.categories) %in% gene.vector,]

  catg_matrix <- matrix(nrow=ncol(IR.categories), ncol=2)
  row.names(catg_matrix) <- colnames(IR.categories)
  colnames(catg_matrix) <- c("Category","No.of.Genes.in.Category")
  catg_matrix[,1]<- colnames(IR.categories)
  
  catg_matrix[,2] <- apply(IR.categories, 2, function(x){length(x[x=="+"])})
  
  columns <- c("No.DEGs.in.category","OR (CI(95%))", "p.value", "p.adj", "DEGs.in.category")
  resultmat <- matrix(ncol=length(columns), nrow=nrow(catg_matrix))
  colnames(resultmat)<- columns
  row.names(resultmat) <- row.names(catg_matrix)
  ORA.test.res <- as.data.frame(cbind(catg_matrix, resultmat))
  
  for (i in 1:nrow(ORA.test.res)){
    category <- ORA.test.res[i,1]
    no.of.genes <- as.numeric(ORA.test.res$No.of.Genes.in.Category[i])
    among.sign_BOTH <- nrow(panel_geneList_sub[panel_geneList_sub[,category]=="+",])
    ORA.test.res$DEGs.in.category[i] <- paste(row.names(panel_geneList_sub[panel_geneList_sub[,category]=="+",]), collapse = "/")
    ORA.test.res$No.DEGs.in.category[i] <- paste0(among.sign_BOTH, " (out of ", length(gene.vector),")")
    
    if (among.sign_BOTH==0){
      ORA.test.res[i,c("OR (CI(95%))", "p.value")] <- NA
    } else {
      fish.res <- perform_fisher(crosstab(among.sign_BOTH, gene.vector,no.of.genes))
      ORA.test.res[i,c("OR (CI(95%))", "p.value")] <- c(paste0(round(as.numeric(fish.res[2]),4)," (",fish.res[1],")"), fish.res[3])
    }
  }
  ORA.test.res$p.adj <- p.adjust(ORA.test.res$p.value, method="bonferroni")
  colnames(ORA.test.res)[3:ncol(ORA.test.res)] <- paste(Name, colnames(ORA.test.res)[3:ncol(ORA.test.res)], sep="_")
  return(ORA.test.res)
}


