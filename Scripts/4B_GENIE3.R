library("tximport")
library("GENIE3")
library("dplyr")


###building transcript names
t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")
#placenta <- read.table("Z:/hhvu/placentaGenes.txt")

###import count files e7.5 =====
sample <- c("e7.5_2.tsv", "e7.5_3.tsv", "e7.5_4.tsv", "e7.5_5.tsv", "e7.5_6.tsv")
new_E7.5Files <- file.path("Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/", sample)

e7.5.kallisto.tsv <- tximport(new_E7.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")

###GENIE3
exprMatr <- as.matrix(e7.5.kallisto.tsv$counts)
colnames(exprMatr) <- c("e7.5_2", "e7.5_3", "e7.5_4", "e7.5_5", "e7.5_6") #e7.5
exprMatr <- t(scale(t(exprMatr)))

regulators <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_TFensGenes.txt", header = F)
regulators <- regulators[regulators$V1 %in% rownames(exprMatr), "V1"]

set.seed(123) # For reproducibility of results
weightMat <- GENIE3(exprMatr, regulators = as.vector(unique(regulators)))

dim(weightMat)
weightMat[1:5,1:5]

linkList <- getLinkList(weightMat)
dim(linkList)

linkList2 <- getLinkList(weightMat, threshold=quantile(linkList$weight, 0.9))
dim(linkList2)

write.table(linkList2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/e7.5/e7.5_0.9.txt", row.names = F, quote = F, sep = "\t")


###import count files e8.5 =====
sample <- c("e8.5_1.tsv", "e8.5_2.tsv", "e8.5_3.tsv", "e8.5_4.tsv", "e8.5_5.tsv", "e8.5_6.tsv")
new_E8.5Files <- file.path("Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/", sample)

e8.5.kallisto.tsv <- tximport(new_E8.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")

###GENIE3
exprMatr <- as.matrix(e8.5.kallisto.tsv$counts)
colnames(exprMatr) <- c("e8.5_1", "e8.5_2", "e8.5_3", "e8.5_4", "e8.5_5", "e8.5_6") #e8.5
exprMatr <- t(scale(t(exprMatr)))

regulators <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_TFensGenes.txt", header = F)
regulators <- regulators[regulators$V1 %in% rownames(exprMatr), "V1"]

set.seed(123) # For reproducibility of results
weightMat <- GENIE3(exprMatr, regulators = as.vector(unique(regulators)))

dim(weightMat)
weightMat[1:5,1:5]

linkList <- getLinkList(weightMat)
dim(linkList)

linkList2 <- getLinkList(weightMat, threshold=quantile(linkList$weight, 0.9))
dim(linkList2)

write.table(linkList2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/e8.5/20220324_e8.5_0.9.txt", row.names = F, quote = F, sep = "\t")


###import count files e9.5 =====
sample <- c("e9.5_1.tsv", "e9.5_2.tsv", "e9.5_3.tsv", "e9.5_4.tsv", "e9.5_5.tsv")
new_e9.5Files <- file.path("Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/", sample)

e9.5.kallisto.tsv <- tximport(new_e9.5Files, type = "kallisto", tx2gene = t2g, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")

###GENIE3
exprMatr <- as.matrix(e9.5.kallisto.tsv$counts)
colnames(exprMatr) <- c("e9.5_1", "e9.5_2", "e9.5_3", "e9.5_4", "e9.5_5") #e9.5
exprMatr <- t(scale(t(exprMatr)))

regulators <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_TFensGenes.txt", header = F)
regulators <- regulators[regulators$V1 %in% rownames(exprMatr), "V1"]

set.seed(123) # For reproducibility of results
weightMat <- GENIE3(exprMatr, regulators = as.vector(unique(regulators)))

dim(weightMat)
weightMat[1:5,1:5]

linkList <- getLinkList(weightMat)
dim(linkList)

linkList2 <- getLinkList(weightMat, threshold=quantile(linkList$weight, 0.9))
dim(linkList2)

write.table(linkList2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/e9.5/20220324_e9.5_0.9.txt", row.names = F, quote = F, sep = "\t")
