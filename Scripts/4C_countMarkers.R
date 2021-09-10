t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")

e7.5specific <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_ensGenes.txt", header = F)
e7.5hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e7.5Group.txt", header = T)
e8.5specific <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_ensGenes.txt", header = F)
e8.5hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e8.5Group.txt", header = T)
e9.5specific <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_ensGenes.txt", header = F)
e9.5hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e9.5Group.txt", header = T)

#trophoblast giant cell differentiation
tgc <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/0_e7.5_TGCdiffMarkers.txt", header = F)
tgc <- dplyr::inner_join(tgc, t2g[,2:3], by = c("V1" = "ext_gene"))
tgc <- dplyr::distinct(tgc)

tgc[tgc$ens_gene %in% e7.5specific$V1,]
tgc[tgc$ens_gene %in% e7.5hier$ens_gene,]

tgc[tgc$ens_gene %in% e8.5specific$V1,]
tgc[tgc$ens_gene %in% e8.5hier$ens_gene,]
dim(tgc[tgc$ens_gene %in% e8.5hier$ens_gene,])

tgc[tgc$ens_gene %in% e9.5specific$V1,]
tgc[tgc$ens_gene %in% e9.5hier$ens_gene,]
dim(tgc[tgc$ens_gene %in% e9.5hier$ens_gene,])


#ectoplacental cone and spongiotrophoblast maintainance
epcSpt <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/0_e7.5_EPC-SGT-markers.txt", header = F)
epcSpt <- dplyr::inner_join(epcSpt, t2g[,2:3], by = c("V1" = "ext_gene"))
epcSpt <- dplyr::distinct(epcSpt)

epcSpt[epcSpt$ens_gene %in% e7.5specific$V1,]
epcSpt[epcSpt$ens_gene %in% e7.5hier$ens_gene,]
dim(epcSpt[epcSpt$ens_gene %in% e7.5hier$ens_gene,])

epcSpt[epcSpt$ens_gene %in% e8.5specific$V1,]
epcSpt[epcSpt$ens_gene %in% e8.5hier$ens_gene,]
dim(epcSpt[epcSpt$ens_gene %in% e8.5hier$ens_gene,])

epcSpt[epcSpt$ens_gene %in% e9.5specific$V1,]
epcSpt[epcSpt$ens_gene %in% e9.5hier$ens_gene,]
dim(epcSpt[epcSpt$ens_gene %in% e9.5hier$ens_gene,])

#Trophoblast stem cell specification and maintenance
tsc <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/0_e7.5_TSCspecificationMarkers.txt", header = F)
tsc <- dplyr::left_join(tsc, t2g[,2:3], by = c("V1" = "ext_gene"))
tsc <- dplyr::distinct(tsc)

tsc[tsc$ens_gene %in% e7.5specific$V1,]
tsc[tsc$ens_gene %in% e7.5hier$ens_gene,]
dim(tsc[tsc$ens_gene %in% e7.5hier$ens_gene,])

tsc[tsc$ens_gene %in% e8.5specific$V1,]
tsc[tsc$ens_gene %in% e8.5hier$ens_gene,]
dim(tsc[tsc$ens_gene %in% e8.5hier$ens_gene,])

tsc[tsc$ens_gene %in% e9.5specific$V1,]
tsc[tsc$ens_gene %in% e9.5hier$ens_gene,]
dim(tsc[tsc$ens_gene %in% e9.5hier$ens_gene,])

#Chorioallantoic attachment
chorioAll <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/0_e8.5_chorioallantoicAttachment.txt", header = F)
chorioAll <- dplyr::left_join(chorioAll, t2g[,2:3], by = c("V1" = "ext_gene"))
chorioAll <- dplyr::distinct(chorioAll)

chorioAll[chorioAll$ens_gene %in% e7.5specific$V1,]
chorioAll[chorioAll$ens_gene %in% e7.5hier$ens_gene,]
dim(chorioAll[chorioAll$ens_gene %in% e7.5hier$ens_gene,])

chorioAll[chorioAll$ens_gene %in% e8.5specific$V1,]
chorioAll[chorioAll$ens_gene %in% e8.5hier$ens_gene,]
dim(chorioAll[chorioAll$ens_gene %in% e8.5hier$ens_gene,])

chorioAll[chorioAll$ens_gene %in% e9.5specific$V1,]
chorioAll[chorioAll$ens_gene %in% e9.5hier$ens_gene,]
dim(chorioAll[chorioAll$ens_gene %in% e9.5hier$ens_gene,])

#Labyrinth branching and vascularization - Syncytiotrophoblast
laby <- read.table("Z:/hhvu/Project1_2/RNA-seq/Documentation/20210701_newSpecificity-codeCheck/0_e9.5_branching-labyrinthMarkers.txt", header = F)
laby <- dplyr::left_join(laby, t2g[,2:3], by = c("V1" = "ext_gene"))
laby <- dplyr::distinct(laby)

laby[laby$ens_gene %in% e7.5specific$V1,]
laby[laby$ens_gene %in% e7.5hier$ens_gene,]
dim(laby[laby$ens_gene %in% e7.5hier$ens_gene,])

laby[laby$ens_gene %in% e8.5specific$V1,]
laby[laby$ens_gene %in% e8.5hier$ens_gene,]
dim(laby[laby$ens_gene %in% e8.5hier$ens_gene,])

laby[laby$ens_gene %in% e9.5specific$V1,]
laby[laby$ens_gene %in% e9.5hier$ens_gene,]
dim(laby[laby$ens_gene %in% e9.5hier$ens_gene,])
