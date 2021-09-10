library("GenomicFeatures")
library("dplyr")

txdb <- makeTxDbFromGFF("Z:/hhvu/Project1_2/RNA-seq/Mus_musculus.GRCm38.98.gtf", format="gtf")
genes <- genes(txdb)
coor <- as.data.frame(genes)

file <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_ensGenes.txt", header = F)

newCoor <- filter(coor, coor$gene_id %in% file$V1)

# get gene region
for (i in seq(from=1, to=dim(newCoor)[1], by=1)){
  if (newCoor$strand[i] == "+") {
    newCoor$end[i] = newCoor$start[i]
    newCoor$start[i] = newCoor$start[i] - 500
  }
  else if (newCoor$strand[i] == "-") {
    newCoor$start[i] = newCoor$end[i]
    newCoor$end[i] = newCoor$end[i] + 500
  }
}
#manipulate region table
newCoor <- cbind(newCoor, a = 0)
newCoor <- newCoor[, c(1,2,3,6,7,5)] #strand info included
newCoor$seqnames <- paste0("chr", newCoor$seqnames)

write.table(newCoor, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_TSS.bed", col.names = T, row.names = F, quote = F, sep = "\t")
#this file's header will need to be removed before running GREAT


#after running GREAT, load the results here to filter and get the significant terms
#bio process terms
great <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_BP_all.tsv", header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_BP_filtered.csv", quote = F, sep = ',', row.names = F)

#single KO terms
great <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_singleKO_all.tsv", header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_singleKO_filtered.csv", quote = F, sep = ',', row.names = F)

