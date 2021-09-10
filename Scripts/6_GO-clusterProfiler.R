library("clusterProfiler")
library("org.Mm.eg.db")

#STRING e7.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e7.5/largestComponent/e7.5_1/e7.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e7.5_1_BP_filtered.rda")
#dim(goTerms)
grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation
#"positive regulation of immune system process"
#"regulation of cell migration"
#"positive regulation of cell population proliferation"
#"epithelial cell migration"
#"regulation of endothelial cell proliferation"
#"smooth muscle cell proliferation"


#STRING e8.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_1/e8.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_1_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation
#[1] "regulation of cell motility"
#[2] "chordate embryonic development"
#[3] "tube morphogenesis"
#[4] "placenta development"
#[5] "embryonic organ development"
#[6] "vasculature development"
#[7] "epithelial cell proliferation"
#[8] "ameboidal-type cell migration"
#[9] "regulation of embryonic development"
#[10] "morphogenesis of a branching structure"


#STRING e8.5_2
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_2/e8.5_2_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_2_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation


#STRING e8.5_3
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e8.5/largestComponent/e8.5_3/e8.5_3_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e8.5_3_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation
#no interesting terms here with clusterProfiler. Originally with GREAT only interesting singleKO terms

#STRING e9.5_1
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_1/e9.5_1_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_1_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation


#STRING e9.5_2
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_2/e9.5_2_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_2_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation

#STRING e9.5_3
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_3/e9.5_3_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_3_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation
#no interesting terms here with clusterProfiler. Originally with GREAT only interesting singleKO terms

#STRING e9.5_4
geneList <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/e9.5/largestComponent/e9.5_4/e9.5_4_TSS.bed", header = F, sep = "\t")
goTerms <- go(geneList$V4)
save(goTerms, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/STRING/e9.5_4_BP_filtered.rda")
#dim(goTerms)
#grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", goTerms$Description, value = T)
#test the following terms with permutation

