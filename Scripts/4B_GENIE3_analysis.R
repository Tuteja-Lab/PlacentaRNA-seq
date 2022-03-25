library("dplyr")

t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")
placenta <- read.table("Z:/hhvu/placentaGenes.txt")
placenta <- inner_join(placenta, t2g[,2:3], by = c("V1" = "ext_gene"))
placenta <- distinct(placenta)


### load networks - all gene analysis ####
for (time in c("e8.5", "e9.5")) {
  coordinates <- read.table(paste0("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/", time, "specific_TSS.bed"), header = F, sep = "\t")
  dir <- paste0("Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/", time, "/")
  for (i in 1:3) {
    subname <- paste0(time, "_", i)
    network <- read.table(paste0(dir, subname, "/", subname, "_nodeTable.csv"), sep = ",", header = T)
    network <- inner_join(network, t2g[,2:3], by = c("name" = "ens_gene"))
    network <- distinct(network)
    sub_coor <- subset(coordinates, coordinates$V4 %in% network$name)
    write.table(sub_coor, paste0(dir, subname, "/", subname, "_TSS.bed"), sep = "\t", quote = F, row.names = F)
  }
}

for (time in c("e7.5", "e8.5", "e9.5")) {
  for (i in 1:3) {
    geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5B_GENIE3/", time, "/", time, "_", i, "/", time, "_", i, "_TSS.bed"), header = F, sep = "\t")
    goTerms <- go(geneList$V4)
    save(goTerms, file = paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/GENIE3/max1000-noSimplify/", time, "_", i, "_BP_filtered.rda"))
    write.table(goTerms, file = paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/GENIE3/max1000-noSimplify/", time, "_", i, "_BP_filtered.txt"),
                quote = F, sep = "\t", row.names = F)
  }
}

### network analysis ####
for (time in c("e9.5")) {
  coordinates <- read.table(paste0("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/", time, "specific_TSS.bed"), header = F, sep = "\t")
  dir <- paste0("Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/", time, "/")
  for (i in 1:3) {
    subname <- paste0(time, "_", i)
    network <- read.table(paste0(dir, subname, "/", subname, "_nodeTable.csv"), sep = ",", header = T)
    network <- inner_join(network, t2g[,2:3], by = c("name" = "ens_gene"))
    network <- distinct(network)
    print(subname)
    
    print(network[which(network$name %in% placenta$ens_gene), "ext_gene"])
    
    topDegree <- subset(network, network$Degree > quantile(network$Degree, 0.9))
    print(topDegree$ext_gene)
    topCloseness <- subset(network, network$ClosenessCentrality > quantile(network$ClosenessCentrality, 0.9)) 
    print(topCloseness$ext_gene)
    topBetweenness <- subset(network, network$BetweennessCentrality > quantile(network$BetweennessCentrality, 0.9)) 
    print(topBetweenness$ext_gene)
    
    all <- intersect(intersect(topDegree$ext_gene, topCloseness$ext_gene), topBetweenness$ext_gene)
    print(all)
    save(all, file=paste0(dir, subname, "/", subname, "_hubGenes.rda"))
    }
}

