# First, we load the libraries we need to carry out the analysis and generate the plots.
# using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    # Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    # sequence count data: removing the noise and preserving large differences.
    # Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

# @@    http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

library("DESeq2")
library("scales")
library("RColorBrewer")
library("gplots")
library(ggplot2)
library("calibrate")
library("biomaRt")
# We assign the name of the conditions to variables. This way we can use the same script for other analyses by doing some minor changes.
options(stringsAsFactors = FALSE)

s1<-"HUSH-batch"

extension<-"r20"

# Reading the metadata file:
meta <- read.delim(file = paste("meta_", s1,".txt", sep=""), header = F)
sampleTable<-data.frame(sampleName=meta$V1, fileName=meta$V2, condition=meta$V3, experiment=meta$V4)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "expressed_instances_r20", design= ~ condition)
colData(ddsHTSeq)$condition <- relevel(colData(ddsHTSeq)$condition, "shC_1")

# To runDESeq, we just need one line:
ddsHTSeq <- DESeq(ddsHTSeq)

resultsNames(ddsHTSeq)


################################################

ctrl<-"shC"
cnd1<-"shM"

dirname<-paste0(cnd1, "_", ctrl)

coef<-paste0("condition_",cnd1,"_vs_",ctrl)
resLFC <- lfcShrink(ddsHTSeq, coef=coef, type="apeglm")
dir.create(file.path(getwd (), dirname), showWarnings = FALSE)

de.filtered = resLFC[!is.na(resLFC$padj), ]



### Tables
# We can already print the table with the results
write.table (resLFC, file = paste (extension,"_",cnd1, "-", ctrl, "_apeglm.txt", sep = ""), quote = F, row.names = T, col.names = T, sep = "\t")

# We can subset out results to include only significantly, differentially expressed genes.
significant.genes=de.filtered[(de.filtered$padj < .05),]

# Se can also use a threshold for a minimum fold change:
up<- significant.genes[(significant.genes$log2FoldChange > 0),]
down<- significant.genes[(significant.genes$log2FoldChange < 0),]

# We write the filtered results to a plain text file. These will be "tab" delimited tables:
write.table(up,file=paste(extension,"_up_", cnd1, "-", ctrl, ".apeglm.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
write.table(down,file=paste(extension,"_down_", cnd1, "-", ctrl, ".apeglm.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")

### checking by log2FC
up2<- significant.genes[(significant.genes$log2FoldChange > 1),]
down2<- significant.genes[(significant.genes$log2FoldChange < (-1)),]

# We write the filtered results to a plain text file. These will be "tab" delimited tables:
write.table(up2,file=paste("log2FC1/", extension, "_", "up_", cnd1, "-", ctrl, ".apeglm.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
write.table(down2,file=paste("log2FC1/", extension, "_", "down_", cnd1, "-", ctrl, ".apeglm.txt", sep=""), quote=F, row.names=T, col.names=T, sep = "\t")
