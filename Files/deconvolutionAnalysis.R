library("linseed", suppressMessages())
t2g <- read.table("Files/t2g.txt", header = T, sep = "\t")
load("Files/e7.5geneLevelTPM.rda")
load("Files/e8.5geneLevelTPM.rda")
load("Files/e9.5geneLevelTPM.rda")
abundance <- dplyr::inner_join(e7.5abundance, e8.5abundance, by = c("ens_gene" = "ens_gene"))
abundance <- dplyr::inner_join(abundance, e9.5abundance, by = c("ens_gene" = "ens_gene"))
rownames(abundance) <- abundance$ens_gene
keep <- colnames(abundance)[2:17]
abundance <- abundance[,keep]
t2g2 <- dplyr::distinct(t2g[,2:3])
t2g2 <- subset(t2g2, t2g2$ens_gene %in% rownames(abundance))
rownames(t2g2) <- t2g2$ens_gene

lo <- LinseedObject$new(abundance, topGenes=5000, annotation = t2g2, geneSymbol = "ens_gene") #use top 5000 most expressed genes
lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(10000) #use 100000 times to calculate significance. To reduce running time, reduce this number, but the results will be different
temp <- lo$spearman
names(lo$genes$pvals) <- rownames(temp) #ensure the gene names between objects matched

lo$significancePlot(0.05)
lo$filterDatasetByPval(0.05)

pdf('Files/svdPlot.pdf', width = 7, height = 7)
lo$svdPlot()
dev.off()

#calculate variance SVD for full data
dataFull <- as.data.frame(lo$exp$full)
dataFull <- dataFull[, grep("norm", colnames(dataFull))]
vars <- svd(dataFull)$d^2
vars <- cumsum(vars / sum(vars))
dfFull <- data.frame(dimension=1:length(vars), variance=vars, filtered=FALSE)
dfFull
#calculate variance SVD for filtered data
dataFiltered <- as.data.frame(lo$exp$filtered)
dataFiltered <- dataFiltered[, grep("norm", colnames(dataFiltered))]
vars <- svd(dataFiltered)$d^2
vars <- cumsum(vars / sum(vars))
dfFiltered <- data.frame(dimension=1:length(vars), variance=vars, filtered=TRUE)
dfFiltered

set.seed(123)
lo$setCellTypeNumber(5)
lo$project("full") # projecting full dataset
lo$project("filtered")
lo$smartSearchCorners(dataset="filtered", error="norm")
lo$deconvolveByEndpoints()
pdf('Files/proportions.pdf', width = 7, height = 7)
print(plotProportions(lo$proportions))
dev.off()

cells <- lo$signatures
lo$selectGenes(100)

distances <- lo$distances
colnames(distances) <- c(paste0("Group", 1:5))
for (i in 1:5) {
  sub <- distances[lo$markers[[i]],]
  keep <- c()
  for (j in 1:nrow(sub)) {
    if (sub[j, i] == min(sub[j,])) {
      keep <- c(keep, rownames(sub)[j])
    }
  }
  markers <- unique(t2g[t2g$ens_gene %in% keep, "ext_gene"])
  print(markers)
  #write.table(markers, paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/cell", i, "markers-minDistance-inTop100.txt"), quote = F, row.names = F)
  write.table(markers, paste0("Files/cell", i, "markers-minDistance-inTop100.txt"), quote = F, row.names = F)
  
}


e7.5specific <- read.table("Files/e7.5specific_ensGenes.txt")
e8.5specific <- read.table("Files/e8.5specific_ensGenes.txt")
e9.5specific <- read.table("Files/e9.5specific_ensGenes.txt")
for (i in 1:5) {
  markers <- read.table(paste0("Files/cell", i, "markers-minDistance-inTop100.txt"), header = T, sep = "\t")
  print(paste0("Cell group ", i))
  print("Number of markers specific to e7.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e7.5specific$V1), "ext_gene"])))
  print("Number of markers specific to e8.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e8.5specific$V1), "ext_gene"])))
  print("Number of markers specific to e9.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e9.5specific$V1), "ext_gene"])))
}


