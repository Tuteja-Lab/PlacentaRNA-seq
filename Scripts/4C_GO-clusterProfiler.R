library("clusterProfiler")
library("org.Mm.eg.db")

#STRING e7.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e7.5/largestComponent/e7.5_1/e7.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e7.5_1_BP_filtered.rda")

#STRING e8.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_1/e8.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_1_BP_filtered.rda")

#STRING e8.5_2
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_2/e8.5_2_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_2_BP_filtered.rda")

#STRING e8.5_3
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_3/e8.5_3_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_3_BP_filtered.rda")

#STRING e9.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_1/e9.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_1_BP_filtered.rda")

#STRING e9.5_2
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_2/e9.5_2_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_2_BP_filtered.rda")

#STRING e9.5_3
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_3/e9.5_3_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_3_BP_filtered.rda")

#STRING e9.5_4
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_4/e9.5_4_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_4_BP_filtered.rda")
