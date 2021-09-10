library("clusterProfiler")
#library("org.Mm.eg.db")
library("ggplot2")
library("doParallel")

dir <- "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/"
timepoint <- "e7.5"
subnet <- "e7.5_1"
load(paste0(dir, subnet, "_BP_filtered.rda"))
original <- goTerms
interestingTerms <- original[grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description),]

geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/", timepoint, "/largestComponent/", subnet, "/", subnet, "_TSS.bed"), header = F, sep = "\t")
allGenes <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/4_specificGenes/", timepoint, "specific_ensGenes.txt"), header = T)

for (k in 1:dim(interestingTerms)[1]) {
  term <- interestingTerms[k,]
  permRes <- perm(geneList$V4, allGenes$x, 10000, term$Description)
  permRes <- as.data.frame(permRes)
  write.table(permRes, paste0(dir, subnet, "_", term$Description, "_qval.txt"), row.names = F, quote = F)
  pdf(paste0(dir, subnet, "_", term$Description, "_hist.pdf"), height = 5, width = 5)
  ggplot(permRes, aes(x = log10(permRes))) + geom_histogram() +
    geom_vline(data = term, aes(xintercept = -log10(qvalue)), color="blue", linetype="dashed", size=1) +
    labs(title = paste0("\'", term$Description, "\'", "\n", "in ", subnet, " network"), x = "-log10(qvalue)", y = "Count")
  dev.off()
}