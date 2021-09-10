library("dplyr")

t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")
placenta <- read.table("Z:/hhvu/placentaGenes.txt")
placenta <- inner_join(placenta, t2g[,2:3], by = c("V1" = "ext_gene"))
placenta <- distinct(placenta)


### load networks - all gene analysis ####
coordinates <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_TSS.bed", header = F, sep = "\t")
nameMap <- read.table("Z:/hhvu/Project1_2/RNA-seq/5A_STRING/e8.5/e8.5_string_mapping.tsv", sep = "\t", quote = "", header = F)
colnames(nameMap) <- c("queryIndex", "queryItem",	"stringId",	"preferredName",	"annotation")
dir <- "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/e8.5/"
subname <- "e8.5_1"

network <- read.table(paste0(dir, subname, "/", subname, "_nodeTable.csv"), sep = ",", header = T)
network <- inner_join(network, t2g[,2:3], by = c("name" = "ens_gene"))
network <- distinct(network)


sub_coor <- subset(coordinates, coordinates$V4 %in% network$name)
#write.table(sub_coor, paste0(dir, subname, "/", subname, "_TSS.bed"), sep = "\t", quote = F, row.names = F)

#use GREAT for ontology analysis
#then use the next lines for GREAT filtering
#bio process terms
great <- read.table(paste0(dir, subname, "/GREAT/", subname, "_BP_all.tsv"), header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, paste0(dir, subname, "/GREAT/", subname, "_BP_filtered.tsv"), quote = F, sep = ',', row.names = F)

#single KO terms
great <- read.table(paste0(dir, subname, "/GREAT/", subname, "_singleKO_all.tsv"), header = F, sep = "\t", quote = "")
great <- subset(great, great$V3 <= 0.05)
great <- subset(great, great$V4 > 2)
great <- great[, -8]
colnames(great) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
write.table(great, paste0(dir, subname, "/GREAT/", subname, "_singleKO_filtered.tsv"), quote = F, sep = ',', row.names = F)


### network analysis ####
network[which(network$name %in% placenta$ens_gene), "ext_gene"]

topDegree <- subset(network, network$Degree > quantile(network$Degree, 0.9))
topDegree$ext_gene
topCloseness <- subset(network, network$ClosenessCentrality > quantile(network$ClosenessCentrality, 0.9)) 
topCloseness$ext_gene
topBetweenness <- subset(network, network$BetweennessCentrality > quantile(network$BetweennessCentrality, 0.9)) 
topBetweenness$ext_gene

all <- intersect(intersect(topDegree$ext_gene, topCloseness$ext_gene), topBetweenness$ext_gene)
all
save(all, file=paste0(dir, subname, "/", subname, "_hubGenes.rda"))

