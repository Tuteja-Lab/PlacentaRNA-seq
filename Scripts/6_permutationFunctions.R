go <- function(subnetworkGenes) {
  GO <- enrichGO(gene = subnetworkGenes, OrgDb=org.Mm.eg.db, ont = "BP", qvalueCutoff = 0.05, maxGSSize=1000, readable = T, keyType = "ENSEMBL")
  goSimplify <- simplify(GO)
  #goSimplify <- GO
  simplify2 <- as.data.frame(goSimplify)
  simplify2$GeneRatio <- sapply(strsplit(simplify2$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  simplify2$BgRatio <- sapply(strsplit(simplify2$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  simplify2$Fold <- as.numeric(simplify2$GeneRatio)/as.numeric(simplify2$BgRatio)
  simplify2 <- simplify2[simplify2$qvalue <= 0.05 & simplify2$Fold >= 2 & simplify2$Count >= 5,]
  
  return(simplify2)
}

goPerm <- function(permGenes) {
  GO <- enrichGO(gene = permGenes, OrgDb=org.Mm.eg.db, ont = "BP", qvalueCutoff = 1, pvalueCutoff = 1, maxGSSize=1000, readable = T, keyType = "ENSEMBL")
  GO <- as.data.frame(GO)

  return(GO)
}

perm <- function(subnetworkGenes, allGenes, numOfPerm, term) {
  result <- c()
  for (i in 1:numOfPerm) {
    #set.seed(123)
    permGenes <- sample(allGenes, length(subnetworkGenes), replace = T)
    goPerm <- goPerm(permGenes)
    goPermTerm <- goPerm[goPerm$Description == term,]
    if (dim(goPermTerm)[1] == 0) {
      result <- c(result, 1)
    } else {
      result <- c(result, goPermTerm$qvalue)
    }
    print(i)
  }
  
  return(result)
}

