library("tidyverse")

### mapping all genes ####
t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")
#load placenta development genes
placenta <- read.table("Z:/hhvu/placentaGenes.txt")
placenta <- inner_join(placenta, t2g[,2:3], by = c("V1" = "ext_gene"))
placenta <- distinct(placenta)

### load networks - all gene analysis ####
coordinates <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_TSS.bed", header = F, sep = "\t")
nameMap <- read.table("Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/e8.5_string_mapping.tsv", sep = "\t", quote = "", header = F)
colnames(nameMap) <- c("queryIndex", "queryItem",	"stringId",	"preferredName",	"annotation")

network <- read.table("Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_2/e8.5_2_nodeTable.csv", sep = ",", header = T)
network <- inner_join(network, nameMap[,c("queryItem", "preferredName")], by = c("name" = "preferredName"))
network <- inner_join(network, t2g[,2:3], by = c("queryItem" = "ens_gene"))
network <- distinct(network)

sub_coor <- subset(coordinates, coordinates$V4 %in% network$queryItem)
write.table(sub_coor, "Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_3/e8.5_3_TSS.bed", sep = "\t", quote = F, row.names = F)


#use GREAT for ontology analysis
#then use the next lines for GREAT filtering
#bio process terms
great <- read.table("Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_3/GREAT/e8.5_3_BP_all.tsv", header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, "Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_3/GREAT/e8.5_3_BP_filtered.csv", quote = F, sep = ',', row.names = F)

#single KO terms
great <- read.table("Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_3/GREAT/e8.5_3_singleKO_all.tsv", header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, "Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_3/GREAT/e8.5_3_singleKO_filtered.csv", quote = F, sep = ',', row.names = F)


#the next part is to analyze networks and get hub genes
network <- subset(network, network$queryItem %in% coordinates$V4)

network[which(network$queryItem %in% placenta$ens_gene), "ext_gene"] #placenta genes

### network analysis ####
topDegree <- subset(network, network$Degree > quantile(network$Degree, 0.9))
topDegree$ext_gene
topCloseness <- subset(network, network$ClosenessCentrality > quantile(network$ClosenessCentrality, 0.9)) 
topCloseness$ext_gene
topBetweenness <- subset(network, network$BetweennessCentrality > quantile(network$BetweennessCentrality, 0.9)) 
topBetweenness$ext_gene

all <- intersect(intersect(topDegree$ext_gene, topCloseness$ext_gene), topBetweenness$ext_gene)
all

save(all, file="Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/largestComponent/e8.5_2/e8.5_2_hubGenes.rda")


