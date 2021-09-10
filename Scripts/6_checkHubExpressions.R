library(dplyr)
library(xlsx)

t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T)
mouseToHuman <- read.table("Z:/hhvu/mouse-humanOrthologs.txt", header = T, sep = "\t", fill = T)

#get expressions all 3 timepoints
load("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/e7.5geneLevelTPM.rda")
e7.5abundance$meanE7.5 <- rowMeans(e7.5abundance[, 2:ncol(e7.5abundance)])

load("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/e8.5geneLevelTPM.rda")
e8.5abundance$meanE8.5 <- rowMeans(e8.5abundance[, 2:ncol(e8.5abundance)])

load("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/e9.5geneLevelTPM.rda")
e9.5abundance$meanE9.5 <- rowMeans(e9.5abundance[, 2:ncol(e9.5abundance)])

tpm <- inner_join(e7.5abundance[,c(1,ncol(e7.5abundance))], e8.5abundance[,c(1,ncol(e8.5abundance))], by = c("ens_gene" = "ens_gene"))
tpm <- inner_join(tpm, e9.5abundance[,c(1,ncol(e9.5abundance))], by = c("ens_gene" = "ens_gene"))
tpm <- inner_join(tpm, distinct(t2g[,2:3]), by = c("ens_gene" = "ens_gene"))

massFunction <- function(dir, time, sub, type) {
  #load hub genes
  if (type == "STRING") {
    load(paste0(dir, time, "/largestComponent/", sub, "/", sub, "_hubGenes.rda"))
  } else {
    load(paste0(dir, time, "/", sub, "/", sub, "_hubGenes.rda"))
  }
  hub <- distinct(t2g[t2g$ext_gene %in% all, 2:3])
  
  #load Ha's RNA-seq data
  exp1 <- subset(tpm, tpm$ens_gene %in% hub$ens_gene)
  
  #load PCE data
  gene_map <- subset(mouseToHuman, mouseToHuman$Gene.stable.ID %in% hub$ens_gene)
  result <- subset(dataset$PlacentaDeciduaBloodData$expressionData,
                   row.names(dataset$PlacentaDeciduaBloodData$expressionData) %in% gene_map$Gene.stable.ID.1)
  result$Gene <- row.names(result)
  result <- dplyr::inner_join(gene_map, result, by = c("Gene.stable.ID.1" = "Gene"))
  keepCol <- c("Gene.stable.ID", "Gene.stable.ID.1", "Gene.name.1", "EVT", "SCT", "VCT")
  result <- result[, colnames(result) %in% keepCol]
  rownames(result) <- result$Gene.name
  exp2 <- left_join(exp1, result, by = c("ens_gene" = "Gene.stable.ID"))
  
  #D0-7 data
  d0_7 <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/avg_tpm_kallisto_Day0_7.txt")
  colnames(d0_7) <- c("transID", "ext_gene", "D0", "D1", "D2",
                      "D3", "D4", "D5", "D6", "D7")
  d0_7 <- d0_7[d0_7$ext_gene %in% hub$ext_gene,]
  l <- list()
  for (i in 1:length(unique(d0_7$ext_gene))) {
    j <- unique(d0_7$ext_gene)[i]
    keep <- c()
    print(j)
    temp <- d0_7[d0_7$ext_gene == j,]
    m <- max(temp[,3:ncol(temp)])
    if (length(which(temp == m)) > 1) {
      temp <- temp[1,]
    } else {
      temp <- temp[which(temp == m, arr.ind = T)[1],]
    }
    keep <- temp[1,1]
    keep <- c(keep, names(which.max(temp[,3:ncol(temp)])))
    keep <- c(keep, temp[,names(which.max(temp[,3:ncol(temp)]))])
    l[[i]] <- keep
  }
  exp2$D0_7 <- NA
  exp2$D0_7[1:length(l)] <- l
  
  #common human cell lines
  cellLines <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/human.all.single-end.log2tpm.txt", header = T)
  cellLines$meanHTR8 <- rowMeans(cellLines[, grep("HTR8", colnames(cellLines))])
  cellLines$meanBeWo <- rowMeans(cellLines[, grep("BeWo", colnames(cellLines))])
  cellLines$meanJEG3 <- rowMeans(cellLines[, grep("JEG3", colnames(cellLines))])
  
  cellLines <- cellLines[cellLines$aGene_Name %in% exp2$Gene.name.1,]
  
  
  exp2 <- left_join(exp2, cellLines[, c(1, grep("mean", colnames(cellLines)))],
                    by = c("Gene.name.1" = "aGene_Name"))
  
  #Okae data
  okae <- read.csv("Z:/hhvu/Project1_2/RNA-seq/Documentation/Okae_et_al_Table_S1.csv", header = T)
  colnames(okae)[1] <- "Gene"
  okae <- okae[, setdiff(colnames(okae), grep("Stroma", colnames(okae), value = T))]
  okae <- okae[okae$Gene %in% exp2$Gene.name.1,]
  
  f <- function(row) {
    keep <- c()
    m <- max(row[,3:ncol(row)])
    if (length(row[which(row == m, arr.ind = T)]) > 1) {
      keep <- unlist(row[which(row == m)[1]])
    } else {
      keep <- c(keep, unlist(row[which(row == m, arr.ind = T)[2]]))
    }
    return(keep)
  }
  maxOkae <- list()
  for (i in 1:nrow(okae)) {
    t <- f(okae[i,])
    maxOkae[[i]] <- t
  }
  okae$maxOkae <- maxOkae
  
  exp2 <- left_join(exp2, okae[,c("Gene", "maxOkae")], by = c("Gene.name.1" = "Gene"))
  
  #HTR8 Plagl1 KD
  htr8 <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/HTR8-Plagl1KD/2_transcriptLevel-kallisto/allSamples-mappedGenes_TPM.csv", header = T, sep = ",")
  htr8$meanHTR8 <- rowMeans(htr8[,c(10,11,12,13,14,15,16,17)])
  htr8 <- htr8[htr8$geneNames %in% exp2$Gene.name.1,]
  t <- htr8 %>% group_by(geneNames) %>% top_n(1, meanHTR8)
  l <- list()
  for (i in 1:nrow(t)) {
    l[[i]] <- c(t$transcripts[i], t$meanHTR8[i])
  }
  t$newMean <- l
  exp2 <- left_join(exp2, t[,c("geneNames", "newMean")],
                    by = c("Gene.name.1" = "geneNames"))
  exp2[is.na(exp2)] <- "NA"
  
  exp2 <- exp2[,-7]
  colnames(exp2) <- c("mouseEnsGene", "meanE7.5TPM", "meanE8.5TPM", "meanE9.5TPM",
                      "mouseGeneName", "humanEnsGene", "PCE.EVT", "PCE.SCT", "PCE.VCT",
                      "D0_7", "meanHTR8", "meanBeWo", "meanJEG3", "maxOkae", "HTR8.Plagl1KO.ctrl")
  
  exp2 <- exp2[,c("mouseEnsGene", "mouseGeneName", "humanEnsGene",
                  "meanE7.5TPM", "meanE8.5TPM", "meanE9.5TPM",
                  "PCE.EVT", "PCE.SCT", "PCE.VCT",
                  "D0_7", "meanHTR8", "meanBeWo", "meanJEG3", "maxOkae", "HTR8.Plagl1KO.ctrl")]
  return(exp2)
}


load("Z:/hhvu/plancetaCellEnrich.Rdata")

#dir <- "Z:/hhvu/Project1_2/RNA-seq/5A_STRING/"
#type <- "STRING"

dir <- "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/"
type <- "GENIE"

time <- "e9.5"

for (i in c(1,2,3)) {
  sub <- paste0(time, "_", i)
  exp2 <- massFunction(dir, time, sub, type)
  
  write.xlsx2(exp2,
              paste0("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/20210901_hubGenes/", "20210909_", sub, type, "hub.xlsx"),
              sheetName = paste0(sub, type), col.names=T, row.names=F)
}



#snRNA-seq