TF_cofactors <- read.table("Z:/hhvu/mm10_TF_coTF_EnsemblID.txt")
t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")
placenta <- read.table("Z:/hhvu/placentaGenes.txt")
placenta <- dplyr::inner_join(placenta, t2g[,2:3], by = c("V1" = "ext_gene"))
placenta <- dplyr::distinct(placenta)

#e8.5 specific genes
e8.5Hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e8.5Group.txt", header = T)
e8.5_e7.5vsE8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e8.5_e7.5vse8.5_DEtransGenes.txt", header = T)
e8.5_e8.5vsE9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e8.5_e8.5vse9.5_DEtransGenes.txt", header = T)

e8.5specificTrans <- intersect(union(e8.5_e7.5vsE8.5$target_id, e8.5_e8.5vsE9.5$target_id), unique(e8.5Hier$transcripts))
length(e8.5specificTrans)
write.table(e8.5specificTrans, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_trans.txt", row.names = F, quote = F)

e8.5specific <- unique(t2g[t2g$target_id %in% e8.5specificTrans, "ens_gene"])
length(e8.5specific)
write.table(e8.5specific, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_ensGenes.txt", row.names = F, quote = F)

e8.5TransGenes <- t2g[t2g$target_id %in% e8.5specificTrans,]
e7.5Mean <- e8.5Hier[e8.5Hier$transcripts %in% e8.5TransGenes$target_id & e8.5Hier$time == "e7.5",]
e8.5TransGenes <- dplyr::left_join(e8.5TransGenes, e7.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e8.5TransGenes)[ncol(e8.5TransGenes)] <- "e7.5MeanScaled"

e8.5Mean <- e8.5Hier[e8.5Hier$transcripts %in% e8.5TransGenes$target_id & e8.5Hier$time == "e8.5",]
e8.5TransGenes <- dplyr::left_join(e8.5TransGenes, e8.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e8.5TransGenes)[ncol(e8.5TransGenes)] <- "e8.5MeanScaled"

e9.5Mean <- e8.5Hier[e8.5Hier$transcripts %in% e8.5TransGenes$target_id & e8.5Hier$time == "e9.5",]
e8.5TransGenes <- dplyr::left_join(e8.5TransGenes, e9.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e8.5TransGenes)[ncol(e8.5TransGenes)] <- "e9.5MeanScaled"
write.table(e8.5TransGenes, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_TransEnsGenes_meanCtsScaled.txt", row.names = F, quote = F)


e8.5TFs <- intersect(e8.5specific, TF_cofactors$V1)
length(e8.5TFs)
write.table(e8.5TFs, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_TFensGenes.txt", row.names = F, quote = F)

e8.5placenta <- intersect(e8.5specific, placenta$ens_gene)
length(e8.5placenta)
write.table(e8.5placenta, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_placentaEnsGenes.txt", row.names = F, quote = F)

#e7.5 specific genes
e7.5Hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e7.5Group.txt", header = T)
e7.5_e7.5vsE9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e7.5_e7.5vse9.5_DEtransGenes.txt", header = T)

e7.5specificTrans <- setdiff(intersect(e7.5_e7.5vsE9.5$target_id, e7.5Hier$transcripts), e8.5specificTrans)
length(e7.5specificTrans)
write.table(e7.5specificTrans, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_trans.txt", row.names = F, quote = F)

e7.5specific <- unique(t2g[t2g$target_id %in% e7.5specificTrans, "ens_gene"])
length(e7.5specific)
write.table(e7.5specific, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_ensGenes.txt", row.names = F, quote = F)

e7.5TransGenes <- t2g[t2g$target_id %in% e7.5specificTrans,]
e7.5Mean <- e7.5Hier[e7.5Hier$transcripts %in% e7.5TransGenes$target_id & e7.5Hier$time == "e7.5",]
e7.5TransGenes <- dplyr::left_join(e7.5TransGenes, e7.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e7.5TransGenes)[ncol(e7.5TransGenes)] <- "e7.5MeanScaled"

e8.5Mean <- e7.5Hier[e7.5Hier$transcripts %in% e7.5TransGenes$target_id & e7.5Hier$time == "e8.5",]
e7.5TransGenes <- dplyr::left_join(e7.5TransGenes, e8.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e7.5TransGenes)[ncol(e7.5TransGenes)] <- "e8.5MeanScaled"

e9.5Mean <- e7.5Hier[e7.5Hier$transcripts %in% e7.5TransGenes$target_id & e7.5Hier$time == "e9.5",]
e7.5TransGenes <- dplyr::left_join(e7.5TransGenes, e9.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e7.5TransGenes)[ncol(e7.5TransGenes)] <- "e9.5MeanScaled"
write.table(e7.5TransGenes, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_TransEnsGenes_meanCtsScaled.txt", row.names = F, quote = F)

e7.5TFs <- intersect(e7.5specific, TF_cofactors$V1)
length(e7.5TFs)
write.table(e7.5TFs, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_TFensGenes.txt", row.names = F, quote = F)

e7.5placenta <- intersect(e7.5specific, placenta$ens_gene)
length(e7.5placenta)
write.table(e7.5placenta, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_placentaEnsGenes.txt", row.names = F, quote = F)

#e9.5 specific genes
e9.5Hier <- read.table("Z:/hhvu/Project1_2/RNA-seq/2B_hierarchicalClustering/e9.5Group.txt", header = T)
e9.5_e7.5vsE9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e9.5_e7.5vse9.5_DEtransGenes.txt", header = T)

e9.5specificTrans <- setdiff(intersect(e9.5_e7.5vsE9.5$target_id, e9.5Hier$transcripts), e8.5specificTrans)
length(e9.5specificTrans)
write.table(e9.5specificTrans, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_trans.txt", row.names = F, quote = F)

e9.5specific <- unique(t2g[t2g$target_id %in% e9.5specificTrans, "ens_gene"])
length(e9.5specific)
write.table(e9.5specific, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_ensGenes.txt", row.names = F, quote = F)

e9.5TransGenes <- t2g[t2g$target_id %in% e9.5specificTrans,]
e7.5Mean <- e9.5Hier[e9.5Hier$transcripts %in% e9.5TransGenes$target_id & e9.5Hier$time == "e7.5",]
e9.5TransGenes <- dplyr::left_join(e9.5TransGenes, e7.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e9.5TransGenes)[ncol(e9.5TransGenes)] <- "e7.5MeanScaled"

e8.5Mean <- e9.5Hier[e9.5Hier$transcripts %in% e9.5TransGenes$target_id & e9.5Hier$time == "e8.5",]
e9.5TransGenes <- dplyr::left_join(e9.5TransGenes, e8.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e9.5TransGenes)[ncol(e9.5TransGenes)] <- "e8.5MeanScaled"

e9.5Mean <- e9.5Hier[e9.5Hier$transcripts %in% e9.5TransGenes$target_id & e9.5Hier$time == "e9.5",]
e9.5TransGenes <- dplyr::left_join(e9.5TransGenes, e9.5Mean[,c("transcripts", "mean_cts_scaled")], by = c("target_id" = "transcripts"))
colnames(e9.5TransGenes)[ncol(e9.5TransGenes)] <- "e9.5MeanScaled"
write.table(e9.5TransGenes, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_TransEnsGenes_meanCtsScaled.txt", row.names = F, quote = F)

e9.5TFs <- intersect(e9.5specific, TF_cofactors$V1)
length(e9.5TFs)
write.table(e9.5TFs, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_TFensGenes.txt", row.names = F, quote = F)

e9.5placenta <- intersect(e9.5specific, placenta$ens_gene)
length(e9.5placenta)
write.table(e9.5placenta, "Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_placentaEnsGenes.txt", row.names = F, quote = F)
