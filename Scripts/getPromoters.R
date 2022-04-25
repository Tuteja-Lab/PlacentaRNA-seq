library("GenomicFeatures")
library("dplyr")

txdb <- makeTxDbFromGFF("Mus_musculus.GRCm38.98.gtf", format="gtf")
genes <- genes(txdb)
coor <- as.data.frame(genes)

file <- someGeneList #input gene names here
file <- read.table('Files/STRING/e7.5/largestComponent/e7.5_1/e7.5_1_nodeTable_KK.csv', sep= ',', header=T)

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

write.table(newCoor, "file.bed", col.names = T, row.names = F, quote = F, sep = "\t")
#this file's header will need to be removed before running GREAT
