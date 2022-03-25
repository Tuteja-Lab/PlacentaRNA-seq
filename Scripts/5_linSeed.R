library("linseed")
library("tximport")

t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/t2g.txt", header = T, sep = "\t")

sample <- c("S57/abundance.tsv", "S58/abundance.tsv", "S59/abundance.tsv", "S60/abundance.tsv", "S61/abundance.tsv",
            "S50/abundance.tsv", "S62/abundance.tsv", "S63/abundance.tsv", "S64/abundance.tsv", "S65/abundance.tsv", "S66/abundance.tsv",
            "S49/abundance.tsv", "S51/abundance.tsv", "S52/abundance.tsv", "S53/abundance.tsv", "S54/abundance.tsv")
dir <- "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/1_kallisto/official/"
files <- file.path(dir, sample)

e7.5Files <- files[1:5]
e7.5.kallisto <- tximport(e7.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
e7.5abundance <- as.data.frame(e7.5.kallisto$abundance)
e7.5abundance <- tibble::rownames_to_column(e7.5abundance, "ens_gene")
colnames(e7.5abundance) <- c("ens_gene", "e7.5_2", "e7.5_3", "e7.5_4", "e7.5_5", "e7.5_6")
save(e7.5abundance, file=paste0(dir, "e7.5geneLevelTPM.rda"))

e8.5Files <- files[6:11]
e8.5.kallisto <- tximport(e8.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
e8.5abundance <- as.data.frame(e8.5.kallisto$abundance)
e8.5abundance <- tibble::rownames_to_column(e8.5abundance, "ens_gene")
colnames(e8.5abundance) <- c("ens_gene", "e8.5_1", "e8.5_2", "e8.5_3", "e8.5_4", "e8.5_5", "e8.5_6")
save(e8.5abundance, file=paste0(dir, "e8.5geneLevelTPM.rda"))

abundance <- dplyr::inner_join(e7.5abundance, e8.5abundance, by = c("ens_gene" = "ens_gene"))

e9.5Files <- files[12:16]
e9.5.kallisto <- tximport(e9.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
e9.5abundance <- as.data.frame(e9.5.kallisto$abundance)
e9.5abundance <- tibble::rownames_to_column(e9.5abundance, "ens_gene")
colnames(e9.5abundance) <- c("ens_gene", "e9.5_1", "e9.5_2", "e9.5_3", "e9.5_4", "e9.5_5")
save(e9.5abundance, file=paste0(dir, "e9.5geneLevelTPM.rda"))

abundance <- dplyr::inner_join(abundance, e9.5abundance, by = c("ens_gene" = "ens_gene"))
rownames(abundance) <- abundance$ens_gene
keep <- colnames(abundance)[2:17]
abundance <- abundance[,keep]

t2g2 <- dplyr::distinct(t2g[,2:3])
t2g2 <- subset(t2g2, t2g2$ens_gene %in% rownames(abundance))
rownames(t2g2) <- t2g2$ens_gene

lo <- LinseedObject$new(abundance, topGenes=5000, annotation = t2g2, geneSymbol = "ens_gene")

lo$calculatePairwiseLinearity()
lo$calculateSpearmanCorrelation()
lo$calculateSignificanceLevel(100)

temp <- lo$spearman
names(lo$genes$pvals) <- rownames(temp)


pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/linseed.significantPlot.pdf", width = 5, height=5)
lo$significancePlot(0.05)
dev.off()

lo$filterDatasetByPval(0.05)

pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/linseed.svdPlot.pdf", width = 5, height=5)
lo$svdPlot()
dev.off()

save(lo, file="/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/lo.rda")

set.seed(123)
lo$setCellTypeNumber(5)
lo$project("full") # projecting full dataset

pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/linseed.projectionPlot.pdf", width = 5, height=5)
lo$projectionPlot(color="filtered")
dev.off()

lo$project("filtered")
lo$smartSearchCorners(dataset="filtered", error="norm")

lo$deconvolveByEndpoints()

pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/linseed.proportionPlot.pdf", width = 5, height=5)
plotProportions(lo$proportions)
dev.off()


cells <- lo$signatures

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

#list of markers
lo$selectGenes(100)

#get markers
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
}

#check if markers are in time-point specific groups
e7.5specific <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/4_specificGenes/e7.5specific_ensGenes.txt")
e8.5specific <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/4_specificGenes/e8.5specific_ensGenes.txt")
e9.5specific <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/4_specificGenes/e9.5specific_ensGenes.txt")

for (i in 1:5) {
  markers <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/7_deconvolution/cell", i, "markers-minDistance-inTop100.txt"), header = T, sep = "\t")
  print(paste0("Cell group ", i))
  print("Number of markers specific to e7.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e7.5specific$V1), "ext_gene"])))
  print("Number of markers specific to e8.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e8.5specific$V1), "ext_gene"])))
  print("Number of markers specific to e9.5:")
  print(length(unique(t2g[t2g$ens_gene %in% intersect(unique(t2g[t2g$ext_gene %in% markers$x, "ens_gene"]), e9.5specific$V1), "ext_gene"])))
}