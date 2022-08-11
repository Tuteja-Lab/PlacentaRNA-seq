# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

id <- read.table("/work/LAS/geetu-lab/hhvu/soncin.et.al.2018/GPL10558-50081.txt", header = T, fill = T, sep = "\t", quote = "")
idmm10 <- read.table("/work/LAS/geetu-lab/hhvu/soncin.et.al.2018/GPL6885-11608.txt", header = T, fill = T, sep = "\t", quote = "")

# load series and platform data from GEO

#HUMAN
gset <- getGEO("GSE100051", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE100051", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#############
ex2 <- as.data.frame(ex)
ex2$id <- rownames(ex2)
ex2 <- inner_join(ex2, id[, c("ID", "Species", "ILMN_Gene")], by = c("id" = "ID"))

tpm2 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/tpmForClustering.txt", header = T)
t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/t2g.txt", header = T)

load(file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/combine-test-expression1.Rdata")
mouseHumanOrthologs<-dataset$GRCH38$mouseHumanOrthologs
hvGenes <- t2g[t2g$target_id %in% rownames(tpm2),]
hvGenes2 <- mouseHumanOrthologs[mouseHumanOrthologs$Gene.stable.ID %in% hvGenes$ens_gene,]

length(intersect(hvGenes2$Human.gene.name, ex2$ILMN_Gene))
intersect(hvGenes2$Human.gene.name, ex2$ILMN_Gene)

length(unique(hvGenes$ext_gene))
length(unique(hvGenes2$Gene.name))

ex2 <- ex2[ex2$ILMN_Gene %in% hvGenes2$Human.gene.name,]

e7.5 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/e7.5specific_ensGenes.txt", header = F)
e7.5 <- inner_join(e7.5, mouseHumanOrthologs, by = c("V1" = "Gene.stable.ID"))
ex2_e7.5 <- ex2[ex2$ILMN_Gene %in% e7.5$Gene.name,]
ex2_e7.5 <- ex2_e7.5[,1:(ncol(ex2_e7.5)-3)]
ex2_e7.5 <- reshape2::melt(ex2_e7.5)

ggplot(ex2_e7.5, aes(x=variable, y=value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))

e9.5 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/e9.5specific_ensGenes.txt", header = F)
e9.5 <- inner_join(e9.5, mouseHumanOrthologs, by = c("V1" = "Gene.stable.ID"))
ex2_e9.5 <- ex2[ex2$ILMN_Gene %in% e9.5$Gene.name,]
ex2_e9.5 <- ex2_e9.5[,1:(ncol(ex2_e9.5)-3)]
ex2_e9.5 <- reshape2::melt(ex2_e9.5)

ggplot(ex2_e9.5, aes(x=variable, y=value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))



####################

#MOUSE
gset <- getGEO("GSE100052", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6885", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE100052", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

#############
ex2 <- as.data.frame(ex)
ex2$id <- rownames(ex2)
ex2 <- inner_join(ex2, idmm10[, c("ID", "Species", "ILMN_Gene")], by = c("id" = "ID"))

tpm2 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/tpmForClustering.txt", header = T)
t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/t2g.txt", header = T)

hvGenes <- t2g[t2g$target_id %in% rownames(tpm2),]

length(intersect(toupper(hvGenes$ext_gene), ex2$ILMN_Gene))

length(unique(hvGenes$ext_gene))

ex2 <- ex2[ex2$ILMN_Gene %in% toupper(hvGenes$ext_gene),]

e7.5 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/e7.5specific_ensGenes.txt", header = F)
e7.5 <- inner_join(e7.5, distinct(t2g[,2:3]), by = c("V1" = "ens_gene"))
ex2_e7.5 <- ex2[ex2$ILMN_Gene %in% toupper(e7.5$ext_gene),]
ex2_e7.5 <- ex2_e7.5[,1:(ncol(ex2_e7.5)-3)]
ex2_e7.5 <- reshape2::melt(ex2_e7.5)

ggplot(ex2_e7.5, aes(x=variable, y=value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))

e9.5 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/e9.5specific_ensGenes.txt", header = F)
e9.5 <- inner_join(e9.5, distinct(t2g[,2:3]), by = c("V1" = "ens_gene"))
ex2_e9.5 <- ex2[ex2$ILMN_Gene %in% toupper(e9.5$ext_gene),]
ex2_e9.5 <- ex2_e9.5[,1:(ncol(ex2_e9.5)-3)]
ex2_e9.5 <- reshape2::melt(ex2_e9.5)

ggplot(ex2_e9.5, aes(x=variable, y=value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90))


