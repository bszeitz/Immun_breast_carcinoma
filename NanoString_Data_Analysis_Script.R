##############################
#
# NanoString data analysis
#
# Written by Beáta Szeitz
#
# Mainly following the guidelines provided in
# https://doi.org/10.1093/bib/bbaa163 and 
# https://github.com/bhattacharya-a-bt/CBCS_normalization/blob/master/CBCS_normalization_tutorial.pdf
# for data quality assessment and normalization
#
##############################


###########################
# Set wd
###########################

setwd("C:/Users/User/PhD/NanoString/Breast Cancer - Immune Panel/R scripts")


###########################
# Load libraries and custom scripts
###########################

.packages = c("ggbiplot","ggplot2","NanoStringQCPro",
              "MASS","ComplexHeatmap","DESeq2",
              "RUVSeq","ggpubr","readxl", 
              "RColorBrewer", "ggrepel", "gridExtra",
              "cowplot")
lapply(.packages, library, character.only=TRUE)
source("Customized_functions.R")


###########################
# Load annotation table
###########################

Annotation <- read.delim("Annotation.txt")
row.names(Annotation) <- Annotation$SampleID


###########################
# Read RCC files and collect flags across the samples

### Read the files one-by-one, and extract count data and QC flags.
### QC flags are collected in the dataframe "Flag.matrix.sample"
###########################

path <- "C:/Users/User/PhD/NanoString/Breast Cancer - Immune Panel/Raw files"
files.RCC <- paste(path, list.files(path = path), sep="/")
list.files(path = path)


# Create tables that will be filled up based on RCC files
rcc = readRcc(files.RCC[1])
Gene.data = rcc$Code_Summary[,1:3]
row.names(Gene.data) <- Gene.data$Name
raw.expression <- as.data.frame(matrix(nrow=nrow(Gene.data), ncol=length(files.RCC)))
row.names(raw.expression) <- Gene.data$Name
Sample.Data <- as.data.frame(matrix(nrow=length(files.RCC), ncol=11))
colnames(Sample.Data) <- c("File", names(rcc[["Sample_Attributes"]]), "imagingQC", "bindingDensityQC", "limitOfDetectionQC", "positiveLinearityQC")

# Import all RCC files
for (i in 1:length(files.RCC)){
  rcc = readRcc(files.RCC[i])
  raw = rcc$Code_Summary
  print("########## Alignment in gene names: ##########")
  print(all(raw$Name == raw.expression$Name))
  raw.expression[,i] = as.numeric(raw$Count)
  colnames(raw.expression)[i] = strsplit(files.RCC[i],
                                         'Raw files/', fixed=T)[[1]][2]
  Sample.Data[i,1] <- colnames(raw.expression)[i]
  Sample.Data[i,2:7] = as.vector(rcc$Sample_Attributes)
  Sample.Data$Image.Qual[i] = image_quality_flag(rcc)
  Sample.Data$Binding.Dens[i] = binding_density_flag(rcc, 0.1, 2.25) #  high and low specified for MAX/FLEX instruments
  Sample.Data$LoD[i] = LoD_flag(rcc)
  Sample.Data$PosLin[i] = pos_ctrl_lin_flag(rcc)
}
colnames(raw.expression) <- lapply(colnames(raw.expression), FUN=function(x){strsplit(x, split="_")[[1]][3]})

# Create flag matrix for samples
Flag.matrix.sample <- as.data.frame(matrix(nrow=length(files.RCC), ncol=6))
row.names(Flag.matrix.sample) <- Sample.Data$SampleID
colnames(Flag.matrix.sample) <- c("Image.Qual", "Binding.Dens", "LoD", "PosLin", "endogenLoD", "housekLoD")
Flag.matrix.sample[1:4] <- Sample.Data[c("Image.Qual","Binding.Dens","LoD","PosLin")]


###########################
# Define housekeeping, control and endogeneous genes
###########################

housekeeping <- Gene.data[Gene.data$CodeClass =="Housekeeping","Name"]
negative.controls <- Gene.data[Gene.data$CodeClass =="Negative","Name"]
endogeneous <- Gene.data[Gene.data$CodeClass =="Endogenous","Name"]
positive.controls <- Gene.data[Gene.data$CodeClass =="Positive","Name"]


###########################
# Create a flag matrix for the housekeeping genes
###########################

Flag.matrix.housek <- as.data.frame(matrix(nrow=length(housekeeping), ncol=2))
row.names(Flag.matrix.housek) <- housekeeping
colnames(Flag.matrix.housek) <- c("below.Lod", "Subtype.Sign")


###########################
# Define LoD, note number of genes below LoD in the flag matrices
###########################

all(row.names(Flag.matrix.sample) == colnames(raw.expression))
neg.raw <- raw.expression[negative.controls,]
lod = colMeans(neg.raw) - apply(neg.raw,2,function(x){sd(x)}) 

Flag.matrix.sample$endogenLoD <- colSums(raw.expression[row.names(Gene.data)[Gene.data$CodeClass == 'Endogenous'],] < lod)
Flag.matrix.sample$housekLoD <- colSums(raw.expression[row.names(Gene.data)[Gene.data$CodeClass == 'Housekeeping'],] < lod)
Flag.matrix.housek$below.Lod <- rowSums(raw.expression[row.names(Gene.data)[Gene.data$CodeClass == 'Housekeeping'],] < lod)


###########################
# Selection of housekeeping genes
###########################

## 1) No counts below the mean count values of the negative control probes

hk.below.Lod <- row.names(Flag.matrix.housek[Flag.matrix.housek$below.Lod !=0,])
hk.below.Lod
# NUBP1 had counts below LoD.

## 2) No diff. expr. between Ctrl and Met in either subtype

# Create expression table only containing housekeeping genes
subtypes <- c("LUMA", "LUMB1", "TNBC")
hk.raw = raw.expression[housekeeping,row.names(Annotation)]
hk.raw <- apply(hk.raw, c(1,2), as.numeric)

# Create a dataframe containing p-values from differential expression analysis
pval <- as.data.frame(matrix(nrow = nrow(hk.raw), ncol=6))
colnames(pval) <- c("LUMA_MetvsCtrl_p", "LUMB1_MetvsCtrl_p", "TNBC_MetvsCtrl_p", 
                    "LUMA_MetvsCtrl_FDR", "LUMB1_MetvsCtrl_FDR", "TNBC_MetvsCtrl_FDR")
row.names(pval) <- housekeeping

for (i in 1:nrow(hk.raw)){
  for (j in 1:length(subtypes)){
    ctrls <- row.names(Annotation[Annotation$SubtypeCategories == paste0("Ctrl ", subtypes[j]),])
    metas <- row.names(Annotation[Annotation$SubtypeCategories == paste0("Met ", subtypes[j]),])
    Patient.Data.Immune.short_sub <- Annotation[c(ctrls, metas),]
    Patient.Data.Immune.short_sub$SubtypeCategories <- factor(Patient.Data.Immune.short_sub$SubtypeCategories, 
                                                              levels=c(paste0("Ctrl ", subtypes[j]), paste0("Met ", subtypes[j])))
    hk.raw_sub <- hk.raw[,c(ctrls, metas)]
    pval[i,j] = coef(summary(glm.nb(as.numeric(hk.raw_sub[i,]) ~ Patient.Data.Immune.short_sub$SubtypeCategories)))[2,4]
  }
}

# FDR correction on the p-values
pval[,4] <- p.adjust(pval[,1], method="fdr")
pval[,5] <- p.adjust(pval[,2], method="fdr")
pval[,6] <- p.adjust(pval[,3], method="fdr")

# Import the diff. expr. results in the housekeeping genes' flag matrix
Flag.matrix.housek[,2] <- apply(pval[,4:6], 1, FUN=function(x){
  if (any(x < 0.001)){ y <- "***"
  } else if (any(x < 0.01)){ y <- "**"
  } else if (any(x < 0.05)){ y <- "*"
  }else { y <- ""}
})
hk.DE <- row.names(Flag.matrix.housek[Flag.matrix.housek$Subtype.Sign !="",])
hk.DE
# CC2D1B was found to be differentially expressed between Met and Ctrl.


## 3) No outlier mean and CV (assessed based on 1.5xIQR rule)

ratio.matrix <- raw.expression[housekeeping,]
for (i in 1:ncol(ratio.matrix)){
  ratio.matrix[,i] <- log2(raw.expression[housekeeping,i] /  mean(raw.expression[,i], na.rm=T))
}
ratios.sd <- apply(X = ratio.matrix, MARGIN = 1, FUN = stats::sd, na.rm = TRUE)
Flag.matrix.housek$Log2SD <- log2(ratios.sd)
Flag.matrix.housek$Log2Mean <- log2(apply(X = raw.expression[housekeeping,], MARGIN = 1, FUN = mean, na.rm = TRUE))
Flag.matrix.housek$CV <- (2^Flag.matrix.housek$Log2SD) / (2^Flag.matrix.housek$Log2Mean)

CV.thresholds <- c( quantile(Flag.matrix.housek$CV)[2] -1.5*IQR(Flag.matrix.housek$CV) , 
                    quantile(Flag.matrix.housek$CV)[4] +1.5*IQR(Flag.matrix.housek$CV)  )

Mean.thresholds <- c( quantile(Flag.matrix.housek$Log2Mean)[2] -1.5*IQR(Flag.matrix.housek$Log2Mean) , 
                      quantile(Flag.matrix.housek$Log2Mean)[4] +1.5*IQR(Flag.matrix.housek$Log2Mean)  )

Flag.matrix.housek$CV.Sign <- ifelse(Flag.matrix.housek$CV> CV.thresholds[2], "*" ,"" )
Flag.matrix.housek$Mean.Sign <- ifelse(Flag.matrix.housek$Log2Mean < Mean.thresholds[2] & Flag.matrix.housek$Log2Mean > Mean.thresholds[1], "" ,"*" )


# Genes with outlier CV:
hk.outlier.CV <- row.names(Flag.matrix.housek[Flag.matrix.housek$CV.Sign !="",])
hk.outlier.CV
# CC2D1B, NUBP1 and ZKSCAN5 was found to have outlier CV.

# Genes with outlier mean:
hk.outlier.mean <- row.names(Flag.matrix.housek[Flag.matrix.housek$Mean.Sign !="",])
hk.outlier.mean
# CC2D1B, FCF1, TUBB and ZKSCAN5 was found to have outlier mean.


## 4) Good correlation with other housekeeping genes

hk.cormat <- cor(t(raw.expression[housekeeping,]), method = "spearman")
Heatmap(cor(t(hk.cormat), method = "spearman"))

hk.badcorr <- c("ABCF1","NOL7","TUBB","EDC3","MRPS5","GUSB")
# The above specified genes have low correlation (<0.50 Spearman's corr. coeff).


# The selected 30 HK-s that were used for normalization:
housekeeping.to.remove <- unique(c(hk.below.Lod, hk.DE, hk.outlier.CV, hk.outlier.mean, hk.badcorr))

housekeeping.selected <- housekeeping[!housekeeping %in% housekeeping.to.remove]
housekeeping.selected


###########################
# Perform normalization
###########################

set.seed(12345)
RUVg.normalization <- RUVSeq_normalization(Expr = raw.expression, 
                                           Clin = Annotation[colnames(raw.expression),],
                                           Genes = Gene.data, k = 1, 
                                           exclude = housekeeping.to.remove)
Expr.norm.RUVg <- as.data.frame(assay(RUVg.normalization$vsd))

Expr.raw <- raw.expression
Expr.raw[Expr.raw==0] <- 1
Expr.raw <- apply(Expr.raw, c(1,2),FUN=log2)


###########################
# Assess normalization
###########################


## Global gene level variation before and after normalization

MeanvsSD <- as.data.frame(matrix(ncol=3, nrow=nrow(Gene.data)))
row.names(MeanvsSD) <- row.names(Gene.data)
colnames(MeanvsSD) <- c("Log2Mean", "LogSD", "CodeClass")
MeanvsSD$CodeClass <- Gene.data$CodeClass

for (i in 1:nrow(MeanvsSD)){
  MeanvsSD[i,1] <- log2(mean(as.numeric(raw.expression[row.names(MeanvsSD)[i],]), na.rm=TRUE))
  MeanvsSD[i,2] <- log2(sd(as.numeric(raw.expression[row.names(MeanvsSD)[i],]), na.rm=TRUE))
}

ggplot(MeanvsSD, aes(x = Log2Mean, y = LogSD)) +
  geom_point(aes(color = factor(CodeClass))) + ggtitle("Before HK selection")

ggplot(MeanvsSD[MeanvsSD$CodeClass =="Housekeeping",], aes(x = Log2Mean, y = LogSD)) +
  geom_point(aes(color = factor(CodeClass))) + geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.95) +
  xlim(c(min(MeanvsSD[MeanvsSD$CodeClass =="Housekeeping","Log2Mean"]),max(MeanvsSD[MeanvsSD$CodeClass =="Housekeeping","Log2Mean"]))) +
  ylim(c(min(MeanvsSD[MeanvsSD$CodeClass =="Housekeeping","LogSD"]),max(MeanvsSD[MeanvsSD$CodeClass =="Housekeeping","LogSD"])))+ 
  ggtitle("All HK-s")

ggplot(MeanvsSD[housekeeping.selected,], aes(x = Log2Mean, y = LogSD)) +
  geom_point(aes(color = factor(CodeClass))) + geom_smooth(method="loess", se=TRUE, fullrange=FALSE, level=0.95)+ ggtitle("Selected HK-s")

ggplot(MeanvsSD[row.names(MeanvsSD) %in% c(positive.controls, negative.controls, housekeeping.selected, endogeneous),], 
       aes(x = Log2Mean, y = LogSD)) +
  geom_point(aes(color = factor(CodeClass)))+ ggtitle("After HK selection")


for (i in 1:nrow(MeanvsSD)){
  MeanvsSD[i,"Before Norm"] <- sd(as.numeric(Expr.raw[row.names(MeanvsSD)[i],]), na.rm=TRUE) / 
    mean(as.numeric(Expr.raw[row.names(MeanvsSD)[i],]), na.rm=TRUE)
  MeanvsSD[i,"After Norm"] <- sd(as.numeric(Expr.norm.RUVg[row.names(MeanvsSD)[i],]), na.rm=TRUE) / 
    mean(as.numeric(Expr.norm.RUVg[row.names(MeanvsSD)[i],]), na.rm=TRUE)
  
}

MeanvsSD.gg <- reshape2::melt(MeanvsSD[,3:ncol(MeanvsSD)])

ggplot(MeanvsSD.gg, aes(x=value, color=variable)) +
  geom_histogram(aes(y=..density..), alpha=0.1, 
                 position="identity")+
  geom_density(alpha=.2)+xlab("CV")


## PCA plots

for (i in 2:ncol(Annotation)){
  gglist <- list()
  gglist[[1]] <- ggbiplot(prcomp(t(Expr.raw)),
                          circle=F, scale = T, 
                          labels=colnames(Expr.raw), 
                          groups = Annotation[colnames(Expr.raw),colnames(Annotation)[i]],
                          ellipse = T,
                          var.axes	=F) + 
    ggtitle(paste0("Before Norm","\nCategorized by: ",colnames(Annotation)[i]))
  gglist[[2]] <- ggbiplot(prcomp(t(Expr.norm.RUVg)),
                          circle=F, scale = T, 
                          labels=colnames(Expr.norm.RUVg), 
                          groups = Annotation[colnames(Expr.norm.RUVg),colnames(Annotation)[i]],
                          ellipse = T,
                          var.axes	=F) + 
    ggtitle(paste0("After Norm","\nCategorized by: ",colnames(Annotation)[i]))
  print(ggarrange(plotlist=gglist))
}


###########################
# Perform differential expression analysis
###########################

set <- RUVg.normalization$set
pData(set)[["SubtypeCategories"]] <- as.factor(pData(set)[["SubtypeCategories"]])

DESeqDataSet.in <- DESeqDataSetFromMatrix(countData = counts(set)[endogeneous,],
                                          colData = pData(set),
                                          design = ~ 0)

SignifLevel <- 0.05

AllResult <- data.frame(Gene=row.names(Expr.norm.RUVg))

AllResult <- extract_results_DEseq2(AllResult, DESeqDataSet.in, ~  SubtypeCategories,  SignifLevel,
                                    contrast = c("SubtypeCategories","Met LUMA","Ctrl LUMA"), "LUMA_Met_vs_Ctrl")
AllResult <- extract_results_DEseq2(AllResult, DESeqDataSet.in, ~  SubtypeCategories,  SignifLevel,
                                    contrast = c("SubtypeCategories","Met LUMB1","Ctrl LUMB1"), "LUMB1_Met_vs_Ctrl")
AllResult <- extract_results_DEseq2(AllResult, DESeqDataSet.in, ~  SubtypeCategories,  SignifLevel,
                                    contrast = c("SubtypeCategories","Met TNBC","Ctrl TNBC"), "TNBC_Met_vs_Ctrl")

row.names(AllResult) <- AllResult$Gene



###########################
# Overrepresentation analysis
###########################

## Load nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.xls (downloaded from the manufacturer's website)

panel_geneList <- as.data.frame(read_xls("C:/Users/User/PhD/NanoString/Breast Cancer - Immune Panel/nCounter_Human_PanCancer_Immune_Profiling_Panel_Gene_List.xls", sheet=3))
colnames(panel_geneList) <- panel_geneList[1,]
panel_geneList <- panel_geneList[-1,]
row.names(panel_geneList) <- panel_geneList[,1]

## Create a matrix indicating which categories the gene is element of

CTantigens <- na.omit(panel_geneList[panel_geneList$`Gene Class`=="CT Antigen","HUGO Name"])
lack.of.info <- AllResult[!AllResult$Gene %in% panel_geneList$`HUGO Name`,"Gene"]

category.list <- sapply(panel_geneList[,"Immune Response Category"], 
                        function(x){na.omit(strsplit(x, split=", ", fixed=T)[[1]])})
names(category.list) <- row.names(panel_geneList)
category.list <- lapply(category.list, function(x){
  y <- gsub("B-cell Functions", "B-Cell Functions",x )
})
categories <- unique(unlist(category.list))
categories <- categories[order(categories)]

IR.categories <- as.data.frame(matrix(nrow=nrow(panel_geneList), ncol=length(categories)+1))
colnames(IR.categories) <- c(categories, "CT Antigen")
row.names(IR.categories) <- row.names(panel_geneList)

for (i in 1:nrow(IR.categories)){
  genecat <- category.list[[row.names(IR.categories)[i]]]
  if (row.names(IR.categories)[i] %in% CTantigens){
    IR.categories[i,"CT Antigen"] <- "+" 
  } else if (row.names(IR.categories)[i] %in% lack.of.info){
    IR.categories[i,] <- "NA" 
  } else if (length(genecat)!=0){
    IR.categories[i,genecat] <- "+"
  }
  
}
IR.categories[is.na(IR.categories)] <- "-"
IR.categories.colnames <- colnames(IR.categories)

# Filter the matrix for the endogenous genes in our dataset
IR.categories <- merge(IR.categories, AllResult, all=T, by="row.names")
IR.categories <- subset(IR.categories, !is.na(IR.categories$Gene))
row.names(IR.categories) <- IR.categories$Row.names
IR.categories <- IR.categories[,IR.categories.colnames]
IR.categories[is.na(IR.categories)] <- "NA"
colnames(IR.categories)[c(6,8,10)]
colnames(IR.categories)[c(6,8,10)] <- c("Chemokine", "Cytokine","Interleukin")

All.Data <- merge(panel_geneList, AllResult, by="row.names", all.y=T)


## Perform overrepresentation analysis

# Define up- and downregulated genes
LUMA_UP = AllResult[!is.na(AllResult$LUMA_Met_vs_Ctrl_padj) & AllResult$LUMA_Met_vs_Ctrl_padj < 0.05 & 
                      AllResult$LUMA_Met_vs_Ctrl_log2FoldChange >0,"Gene"]
LUMA_DOWN = AllResult[!is.na(AllResult$LUMA_Met_vs_Ctrl_padj) & AllResult$LUMA_Met_vs_Ctrl_padj < 0.05 & 
                        AllResult$LUMA_Met_vs_Ctrl_log2FoldChange <0,"Gene"]
LUMB1_UP = AllResult[!is.na(AllResult$LUMB1_Met_vs_Ctrl_padj) & AllResult$LUMB1_Met_vs_Ctrl_padj < 0.05 & 
                       AllResult$LUMB1_Met_vs_Ctrl_log2FoldChange >0,"Gene"]
LUMB1_DOWN = AllResult[!is.na(AllResult$LUMB1_Met_vs_Ctrl_padj) & AllResult$LUMB1_Met_vs_Ctrl_padj < 0.05 & 
                         AllResult$LUMB1_Met_vs_Ctrl_log2FoldChange <0,"Gene"]
TNBC_UP = AllResult[!is.na(AllResult$TNBC_Met_vs_Ctrl_padj) & AllResult$TNBC_Met_vs_Ctrl_padj < 0.05 & 
                      AllResult$TNBC_Met_vs_Ctrl_log2FoldChange >0,"Gene"]
TNBC_DOWN = AllResult[!is.na(AllResult$TNBC_Met_vs_Ctrl_padj) & AllResult$TNBC_Met_vs_Ctrl_padj < 0.05 & 
                        AllResult$TNBC_Met_vs_Ctrl_log2FoldChange <0,"Gene"]

# Perform the analysis
try.ORA <- ORA_test(LUMA_UP, LUMA_DOWN, "LUMA", IR.categories)
Ncol <- ncol(try.ORA)

ORA.summary <- cbind(ORA_test(LUMA_UP, LUMA_DOWN, "LUMA", IR.categories), 
                     ORA_test(LUMB1_UP, LUMB1_DOWN, "LUMB1", IR.categories)[,c(3:Ncol)],
                     ORA_test(TNBC_UP, TNBC_DOWN, "TNBC", IR.categories)[,c(3:Ncol)])


###########################
# Volcano plot
###########################

## LUMA

de <- AllResult[,c("LUMA_Met_vs_Ctrl_log2FoldChange", "LUMA_Met_vs_Ctrl_padj")]
IRCats <- c("Chemokine", "Complement", "CT Antigen", "Cytokine", "Interleukin")

de <- de[,c(grep("_log2FoldChange", colnames(de)), grep("_padj", colnames(de)))]
colnames(de) <- c("log2FoldChange", "p.adj")
de$p.adj <- ifelse(is.na(de$p.adj), 1, de$p.adj)
de$diffexpressed <- ifelse(de$p.adj < 0.05, "yes", "no")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
DES <- unique(de$delabel)
DES <- DES[DES!=""]
IR.categories.sub <- IR.categories[DES,IRCats]
for (i in 1:ncol(IR.categories.sub)){
  IR.categories.sub[,i] <- gsub("+", colnames(IR.categories.sub)[i], IR.categories.sub[,i], fixed = T)
}
IR.categories.sub[IR.categories.sub=="-"] <- ""
IR.categories.sub$Category <- apply(IR.categories.sub, 1, function(x){
  y <- x[x!=""]
  y <- paste(unique(y), collapse=" + ")})
summary(factor(IR.categories.sub$Category))
de$Category <- ""
for (i in 1:nrow(de)){
  if (de$delabel[i] != ""){
    de$Category[i] <- IR.categories.sub[de$delabel[i],"Category"]
    de$Category[i] <- ifelse(de$Category[i] =="NA", "",de$Category[i])
  } else {
    de$Category[i] <- "Not Significant"
  }
}
de$Category <- ifelse(de$Category =="", "Other", de$Category)

de$Category <- factor(de$Category, levels = c("Chemokine", "Chemokine + Complement", "Chemokine + Cytokine", "Complement",
                                              "CT Antigen", "Cytokine", "Cytokine + Interleukin", "Interleukin", "Other",
                                              "Not Significant"))

title <- "LUMA"

luma <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(p.adj), col=Category, label=delabel)) + 
  geom_point() + 
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text_repel(max.overlaps = 100, seed=1111, mapping = aes(segment.alpha	=0.5)) + ggtitle(title)+ 
  scale_color_manual(values=c("#1F78B4", "#B2DF8A",  "#FB9A99","#33A02C", "#E31A1C", "#FF7F00",
                              "#FDBF6F", "#6A3D9A", "#A6CEE3" ,"darkgrey"))
luma


## LUMB1

de <- AllResult[,c("LUMB1_Met_vs_Ctrl_log2FoldChange", "LUMB1_Met_vs_Ctrl_padj")]
IRCats <- c("CT Antigen", "Leukocyte Functions", "Interleukin")

de <- de[,c(grep("_log2FoldChange", colnames(de)), grep("_padj", colnames(de)))]
colnames(de) <- c("log2FoldChange", "p.adj")
de$p.adj <- ifelse(is.na(de$p.adj), 1, de$p.adj)
de$diffexpressed <- ifelse(de$p.adj < 0.05, "yes", "no")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
DES <- unique(de$delabel)
DES <- DES[DES!=""]
IR.categories.sub <- IR.categories[DES,IRCats]
for (i in 1:ncol(IR.categories.sub)){
  IR.categories.sub[,i] <- gsub("+", colnames(IR.categories.sub)[i], IR.categories.sub[,i], fixed = T)
}
IR.categories.sub[IR.categories.sub=="-"] <- ""
IR.categories.sub$Category <- apply(IR.categories.sub, 1, function(x){
  y <- x[x!=""]
  y <- paste(unique(y), collapse=" + ")})
summary(factor(IR.categories.sub$Category))
de$Category <- ""
for (i in 1:nrow(de)){
  if (de$delabel[i] != ""){
    de$Category[i] <- IR.categories.sub[de$delabel[i],"Category"]
    de$Category[i] <- ifelse(de$Category[i] =="NA", "-",de$Category[i])
  } else {
    de$Category[i] <- "Not Significant"
  }
}
de$Category <- ifelse(de$Category =="", "Other", de$Category)
de$Category <- factor(de$Category, levels = c("CT Antigen", "Interleukin", "Leukocyte Functions","Other", "Not Significant"))

title <- "LUMB1"

lumb1 <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(p.adj), col=Category, label=delabel)) + 
  geom_point() + 
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text_repel(max.overlaps = 100, seed=1111, mapping = aes(segment.alpha	=0.5)) + ggtitle(title)+ 
  scale_color_manual(values=c("#E31A1C", "#6A3D9A","#B15928","#A6CEE3", "darkgrey"))
lumb1


## TNBC

de <- AllResult[,c("TNBC_Met_vs_Ctrl_log2FoldChange", "TNBC_Met_vs_Ctrl_padj")]
IRCats <- c("CT Antigen", "Macrophage Functions")

de <- de[,c(grep("_log2FoldChange", colnames(de)), grep("_padj", colnames(de)))]
colnames(de) <- c("log2FoldChange", "p.adj")
de$p.adj <- ifelse(is.na(de$p.adj), 1, de$p.adj)
de$diffexpressed <- ifelse(de$p.adj < 0.05, "yes", "no")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
de$delabel <- ifelse(de$p.adj < 0.05, row.names(de), "")
DES <- unique(de$delabel)
DES <- DES[DES!=""]
IR.categories.sub <- IR.categories[DES,IRCats]
for (i in 1:ncol(IR.categories.sub)){
  IR.categories.sub[,i] <- gsub("+", colnames(IR.categories.sub)[i], IR.categories.sub[,i], fixed = T)
}
IR.categories.sub[IR.categories.sub=="-"] <- ""
IR.categories.sub$Category <- apply(IR.categories.sub, 1, function(x){
  y <- x[x!=""]
  y <- paste(unique(y), collapse=" + ")})
summary(factor(IR.categories.sub$Category))
de$Category <- ""
for (i in 1:nrow(de)){
  if (de$delabel[i] != ""){
    de$Category[i] <- IR.categories.sub[de$delabel[i],"Category"]
    de$Category[i] <- ifelse(de$Category[i] =="NA", "-",de$Category[i])
  } else {
    de$Category[i] <- "Not Significant"
  }
}
de$Category <- ifelse(de$Category =="", "Other", de$Category)
de$Category <- factor(de$Category, levels = c("CT Antigen", "Macrophage Functions","Other",  "Not Significant"))

title <- "TNBC"

tnbc <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(p.adj), col=Category, label=delabel)) + 
  geom_point() + 
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text_repel(max.overlaps = 100, seed=1111, mapping = aes(segment.alpha	=0.5)) + ggtitle(title)+ 
  scale_color_manual(values=c("#E31A1C", "#CAB2D6","#A6CEE3",  "darkgrey"))
tnbc


top_row <- plot_grid(luma,
                     ncol = 1,
                     labels = c('A'),
                     label_fontfamily = '',
                     label_fontface = 'bold',
                     label_size = 16,
                     align = 'h',
                     rel_widths = c(1))

bottom_row <- plot_grid(lumb1, tnbc,
                        ncol = 2,
                        labels = c('B', 'C'),
                        label_fontfamily = '',
                        label_fontface = 'bold',
                        label_size = 16,
                        align = 'h',
                        rel_widths = c(1, 1))

plot_grid(top_row, bottom_row, ncol = 1,
          rel_heights = c(1.25,0.75))




###########################
# Unsupervised clustering
###########################

Expr.norm.RUVg.scaled <- t(scale(t(Expr.norm.RUVg[endogeneous,row.names(Annotation)])))

#set.seed(1234)
#ht.CE <- Heatmap(Expr.norm.RUVg.scaled, 
#                 row_km = 4, row_km_repeats = 1, clustering_distance_rows = "euclidean",clustering_method_rows = "average")
#ht.CE <- draw(ht.CE)
#c1 <- t(t(row.names(Expr.norm.RUVg.scaled[row_order(ht.CE)[[1]],])))
#c2 <- t(t(row.names(Expr.norm.RUVg.scaled[row_order(ht.CE)[[2]],])))
#c3 <- t(t(row.names(Expr.norm.RUVg.scaled[row_order(ht.CE)[[3]],])))
#c4 <- t(t(row.names(Expr.norm.RUVg.scaled[row_order(ht.CE)[[4]],])))
#save(c1, c2, c3, c4, file="Gene_Clusters.RData")
load("Gene_Clusters.RData") # to make cluster results fixed

gene.order <- c(c1, c2, c3, c4)
clusters.assign <- c(rep("c1", length(c1)), 
                     rep("c2", length(c2)),
                     rep("c3", length(c3)),
                     rep("c4", length(c4)))
names(clusters.assign) <- gene.order

Expr.norm.RUVg.scaled <- Expr.norm.RUVg.scaled[names(clusters.assign),row.names(Annotation)]


row.ha <- rowAnnotation(`Log2\nExpression` = log2(AllResult[row.names(Expr.norm.RUVg.scaled),"LUMA_Met_vs_Ctrl_baseMean"]),
                             col = list(`Log2\nExpression` = circlize::colorRamp2(c(2, 5, 7, 9, 15), c("darkgreen","lightgreen", "yellow", "orange", "brown"))))

column.ha <- HeatmapAnnotation(df=Annotation[,-1],
                               which="col",
                               col=list(
                                 "SubtypeCategories"=c( 
                                   "Met LUMA" = "#33A02C",
                                   "Ctrl LUMA" = "#B2DF8A",
                                   "Met LUMB1" = "#6A3D9A",
                                   "Ctrl LUMB1" = "#CAB2D6",
                                   "Met TNBC" = "#FF7F00", 
                                   "Ctrl TNBC" = "#FDBF6F"
                                 ),
                                 "DMFS"=c(
                                   "[3,30]" = "#deb3b5",
                                   "(30,81]" = "#c17378",
                                   "(81,146]" = "#a3464e",
                                   "(146,198]" = "#6a1f26" 
                                 ),
                                 "strTIL"=c(
                                   "1-3" = "#ffdc74",
                                   "5" = "#ffbf00",
                                   "10-20" = "#c09515"
                                 ),
                                 "N"=c(
                                   "0" = "#C6DBEF",
                                   "1" = "#9ECAE1",
                                   "2-3" = "#2171B5"
                                 ),
                                 "iTIL"=c(
                                   "1" = "#ffdc74",
                                   "2-5" = "#c09515"
                                 ),
                                 "pTIL"=c(
                                   "1" = "#ffdc74",
                                   "2-3" = "#c09515"
                                 ),
                                 "Ki67"=c(
                                   "0-10" = "#beeaac",
                                   "11-29" = "#7ecc67",
                                   "30-80" = "#61934e"
                                 )),
                               gap = unit(c(2.5,2,1,1,1,1), "mm"),show_legend =FALSE
)
col.fun = circlize::colorRamp2(c(-2, 0, 2), c("#0000FF", "#FFFFFF", "#FF0000"))


ht <- Heatmap(Expr.norm.RUVg.scaled, top_annotation = column.ha, 
              name = "Unsupervised Clustering",
              clustering_distance_columns = "euclidean", 
              clustering_distance_rows = "euclidean",
              clustering_method_columns = "average", 
              clustering_method_rows = "average", cluster_column_slices = T,
              col=col.fun,
              column_split = factor(Annotation$SubtypeCategories, 
                                    levels=c("Ctrl LUMA", "Met LUMA", 
                                             "Ctrl LUMB1", "Met LUMB1", 
                                             "Ctrl TNBC", "Met TNBC")),
              show_row_names = F,
              show_row_dend = F,
              row_split = clusters.assign,
              right_annotation = row.ha,
              width = unit(15, "cm"),
              gap = unit(0.25, "cm"),
              column_title_gp = gpar(fontsize = 12),
              row_title_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 8), show_heatmap_legend = F)
ht

#svg("C:/Users/User/PhD/NanoString/Breast Cancer - Immune Panel/Figures/Clustering_Heatmap_vsJan2022.svg", width = 15, height = 10)
#ht
#dev.off()


###########################
# Overrepresentation analysis for gene clusters
###########################

try.ORA <- ORA_test_onlyBOTH(c1, "C1", IR.categories)
Ncol <- ncol(try.ORA)

ORA.summary.clusters <- cbind(ORA_test_onlyBOTH(c1, "C1", IR.categories), 
                              ORA_test_onlyBOTH(c2, "C2", IR.categories)[,c(3:Ncol)],
                              ORA_test_onlyBOTH(c3, "c3", IR.categories)[,c(3:Ncol)], 
                              ORA_test_onlyBOTH(c4, "c4", IR.categories)[,c(3:Ncol)])


###########################
# Export results
###########################

all(ORA.summary[,1] == ORA.summary.clusters[,1])
all(ORA.summary[,2] == ORA.summary.clusters[,2])

ORA.all.results <- cbind(ORA.summary, ORA.summary.clusters[,-c(1:2)])
All.Data.Clusters <- cbind(All.Data, names(clusters.assign[All.Data$Row.names]),clusters.assign[All.Data$Row.names])
colnames(All.Data.Clusters)[31:32] <- c("Gene (cluster)","Cluster")

Expr.norm.RUVg.to.export <- cbind(row.names(Expr.norm.RUVg), Expr.norm.RUVg)
colnames(Expr.norm.RUVg.to.export)[1] <- "Gene"
Expr.norm.RUVg.to.export <- Expr.norm.RUVg.to.export[order(Expr.norm.RUVg.to.export$Gene),]

#write.table(ORA.all.results, "ORA_results.txt", sep="\t", quote = F, row.names = F, na = "")
#write.table(All.Data.Clusters, "Gene_Annotations_DEresults_Clusters.txt", sep="\t", quote = F, row.names = F, na = "")
#write.table(Expr.norm.RUVg.to.export, "Normalized_Expr.txt", sep="\t", quote = F, row.names = F, na = "")


