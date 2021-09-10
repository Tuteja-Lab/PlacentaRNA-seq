library(igraph)
library(xlsx)

t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T)
placenta <- read.table("Z:/hhvu/placentaGenes.txt", header = F)
important <- c("Mmp9", "Tlr2", "Cdh2","Mapk14","Akt1","Mapk1","Hsp90aa1","Creb1",
               "Trip12","Ccnb1","Ccna2","Wnt5a","Fn1","Tgfb1","Vegfa","Egfr","Igf2",
               "Cdh1","Igf1","Col1a1","Spp1","Timp1","Acta2","Gcgr","Gna12","Rhoj",
               "Hmox1","Nr2f2","Ctbp2","Prdm1","Satb1","Id3","Twist1","Tbx4",
               "Tead1","Ss18","Tcf7l1","Bcl9l","Dnmt1","Setd2","Erf","Mycn","Smarcc1",
               "Esx1","Bcl9l","Cebpb","Ncoa3","Tead1","Erf","Cebpa","Tgfb1","Ovol2","Ccnd1",
               "Rb1","Etv3","Elf1","Cdk5","Foxo3","Prnp","Senp1","Jun","Vegfa","Hes1",
               "Cdk2","Zfpm1","Plagl1")


#in STRING directory, do (change timepoint accordingly)
#for i in 1 2 4; do sed 's/,/\t/g' e9.5/largestComponent/e9.5_"$i"/e9.5_"$i"_edgeTable.csv | cut -f 11 | sed 's/\"//g' | sed 's/ (interacts with) /\t/g' | grep -v "name" > e9.5/largestComponent/e9.5_"$i"/e9.5_"$i"_edgeTable_igraph.txt; done
#in GENIE3 directory, do (change timepoint accordingly)
# for i in 1 2 3; do sed 's/,/\t/g' e9.5/e9.5_"$i"/e9.5_"$i"_edgeTable.csv | cut -f 4 | sed 's/\"//g' | sed 's/ (interacts with) /\t/g' | grep -v "name" > e9.5/e9.5_"$i"/e9.5_"$i"_edgeTable_igraph.txt; done

###STRING
dir <- "Z:/hhvu/Project1_2/RNA-seq/5A_STRING/"
time <- "e9.5"
type <- "STRING"

for (i in c(1,2,4)) {
  sub <- paste0(time, "_", i)
  
  file <- paste0(dir, time, "/largestComponent/", sub, "/", sub, "_edgeTable_igraph.txt")
  edge <- read.table(file, header = F, sep = "\t")
  load(paste0(dir, time, "/largestComponent/", sub, "/", sub, "_hubGenes.rda"))
  hub <- all
  graph <- graph_from_data_frame(edge, directed = F)
  subImportant <- unique(placenta[c(which(placenta$V1 %in% edge$V1), which(placenta$V1 %in% edge$V2)), "V1"])
  subImportant <- unique(c(subImportant, intersect(important, edge$V1), intersect(important, edge$V2)))
  
  df <- data.frame(matrix(nrow = length(hub), ncol = length(subImportant)))
  for (j in 1:length(hub)) {
    path <- shortest_paths(graph = graph, from = as.numeric(V(graph)[hub[j]]),
                           to = as.numeric(V(graph)[subImportant]),
                           output = "vpath")
    temp <- as.data.frame(t(sapply(path$vpath, as_ids)))
    if (dim(temp)[1] == 1) {
      df[j,] <- temp
    } else {
      for (k in 1:dim(temp)[1]) {
        df[j,k] <- toString(temp[k,])
      }
    }
  }
  rownames(df) <- hub
  colnames(df) <- subImportant
  
  write.xlsx2(df,
              paste0("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/20210901_hubGenes/", "20210907_", sub, type, "hub-shortestPaths.xlsx"),
              sheetName = paste0(sub, type), col.names=T, row.names=T)
  
  neigh <- c()
  for (k in hub) {
    neigh <- c(neigh, as_ids(neighbors(graph, k)))
  }
  neigh <- unique(neigh)
  df2 <- data.frame(matrix(nrow = length(neigh), ncol = 3))
  df2$X1 <- neigh
  for (l in 1:length(neigh)) {
    df2[l, "X2"] <- toString(intersect(as_ids(neighbors(graph, neigh[l])), hub))
  }
  df2$X3 <- paste0(sub, type)
  colnames(df2) <- c("Gene", "Hubs", "Network")
  
  write.xlsx2(df2,
              paste0("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/20210901_hubGenes/", "20210910_", sub, type, "geneToHubs.xlsx"),
              sheetName = paste0(sub, type), col.names=T, row.names=T)
  
}

###GENIE
dir <- "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/"
time <- "e9.5"
type <- "GENIE"

for (i in c(1,2,3)) {
  sub <- paste0(time, "_", i)
  
  file <- paste0(dir, time, "/", sub, "/", sub, "_edgeTable_igraph.txt")
  edge <- read.table(file, header = F, sep = "\t")
  edge <- dplyr::inner_join(edge, dplyr::distinct(t2g[,2:3]), by = c("V1" = "ens_gene"))
  colnames(edge) <- c("gene1", "gene2", "V1")
  edge <- dplyr::inner_join(edge, dplyr::distinct(t2g[,2:3]), by = c("gene2" = "ens_gene"))
  colnames(edge) <- c("gene1", "gene2", "V1", "V2")
  edge <- edge[,c("V1", "V2")]
  load(paste0(dir, time, "/", sub, "/", sub, "_hubGenes.rda"))
  hub <- all
  graph <- graph_from_data_frame(edge, directed = F)
  subImportant <- unique(placenta[c(which(placenta$V1 %in% edge$V1), which(placenta$V1 %in% edge$V2)), "V1"])
  subImportant <- unique(c(subImportant, intersect(important, edge$V1), intersect(important, edge$V2)))
  
  df <- data.frame(matrix(nrow = length(hub), ncol = length(subImportant)))
  for (j in 1:length(hub)) {
    path <- shortest_paths(graph = graph, from = as.numeric(V(graph)[hub[j]]),
                           to = as.numeric(V(graph)[subImportant]),
                           output = "vpath")
    temp <- as.data.frame(t(sapply(path$vpath, as_ids)))
    if (dim(temp)[1] == 1) {
      df[j,] <- temp
    } else {
      for (k in 1:dim(temp)[1]) {
        df[j,k] <- toString(temp[k,])
      }
    }
  }
  rownames(df) <- hub
  colnames(df) <- subImportant
  
  write.xlsx2(df,
              paste0("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/20210901_hubGenes/", "20210907_", sub, type, "hub-shortestPaths.xlsx"),
              sheetName = paste0(sub, type), col.names=T, row.names=T)
  
  neigh <- c()
  for (k in hub) {
    neigh <- c(neigh, as_ids(neighbors(graph, k)))
  }
  neigh <- unique(neigh)
  df2 <- data.frame(matrix(nrow = length(neigh), ncol = 3))
  df2$X1 <- neigh
  for (l in 1:length(neigh)) {
    df2[l, "X2"] <- toString(intersect(as_ids(neighbors(graph, neigh[l])), hub))
  }
  df2$X3 <- paste0(sub, type)
  colnames(df2) <- c("Gene", "Hubs", "Network")
  
  write.xlsx2(df2,
              paste0("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/20210901_hubGenes/", "20210910_", sub, type, "geneToHubs.xlsx"),
              sheetName = paste0(sub, type), col.names=T, row.names=T)
  
}

path <- shortest_paths(graph = graph, from = as.numeric(V(graph)["Hnrnpk"]),
                       to = as.numeric(V(graph)[hub]),
                       output = "vpath")