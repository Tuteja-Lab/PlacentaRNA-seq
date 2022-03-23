library("dplyr")
library("tidyverse")
library("biomaRt")
library("ggplot2")
library("dendextend")
library("matrixStats")


###building transcript names =====
t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/t2g.txt", header = T, sep = "\t")
t2g <- t2g[order(t2g$target_id),]

#load protein-coding transcripts
coding <- read.table("/work/LAS/geetu-lab/hhvu/Mus_musculus_grcm38_coding_transcripts.txt", header = F)

###filter based on raw est. counts and remove outliers =====
est_counts <- read.csv("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/1_kallisto/official/estCounts_allSamples.tsv", header = TRUE, sep = "\t", stringsAsFactors=FALSE)

#reformat count file
est_counts <- est_counts[order(est_counts$target_id),]
rownames(est_counts) <- est_counts$target_id #make transcript names to be row names
est_counts <- est_counts[, -c(1, 20, 21)] #exclude non-count columns

#remove outliers
est_counts <- est_counts[, -which(names(est_counts) %in% c("E7.5_1", "E9.5_6"))]
E7.5_mean <- rowMeans(est_counts[,c("E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
E8.5_mean <- rowMeans(est_counts[,c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
E9.5_mean <- rowMeans(est_counts[,c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5")])
mean_cts <- data.frame(matrix(ncol = 3, nrow = nrow(est_counts)))
row.names(mean_cts) <- rownames(est_counts)
colnames(mean_cts) <- c("E7.5_mean", "E8.5_mean", "E9.5_mean")
mean_cts$E7.5_mean <- E7.5_mean
mean_cts$E8.5_mean <- E8.5_mean
mean_cts$E9.5_mean <- E9.5_mean

#filter low count transcripts
basic_filter <- function (row, min_reads = 20, min_prop = 1/3) {
  mean(row >= min_reads) >= min_prop
}
keep <- apply(mean_cts, 1, basic_filter)


###load TPM file in to do clustering =====
tpm <- read.csv("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/1_kallisto/official/TPM_allSamples.tsv", header = TRUE, sep = "\t", stringsAsFactors=FALSE)
tpm <- tpm[order(tpm$target_id),]
row.names(tpm) <- tpm$target_id
tpm <- tpm[,-c(1, 20, 21)]
tpm <- tpm[, -which(names(tpm) %in% c("E7.5_1", "E9.5_6"))]
tpm[] <- lapply(tpm, function(x) as.numeric(as.character(x)))
dim(tpm)[1] #number of transcripts originally

#filter out the low count transcripts according to the est_counts filtering above
tpm <- tpm[keep, ]
dim(tpm)[1] #number of transcripts after basic filtering

E7.5_mean <- rowMeans(tpm[,c("E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
E8.5_mean <- rowMeans(tpm[,c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
E9.5_mean <- rowMeans(tpm[,c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5")])
mean_tpm <- data.frame(matrix(ncol = 3, nrow = nrow(tpm)))
row.names(mean_tpm) <- rownames(tpm)
colnames(mean_tpm) <- c("E7.5_mean", "E8.5_mean", "E9.5_mean")
mean_tpm$E7.5_mean <- E7.5_mean
mean_tpm$E8.5_mean <- E8.5_mean
mean_tpm$E9.5_mean <- E9.5_mean
mean_tpm$var <- rowVars(as.matrix(mean_tpm)) #calculate the variance
mean_tpm_keep <- subset(mean_tpm, mean_tpm$var > quantile(mean_tpm$var, 0.25)) #keep the top 75% most variable transcripts

tpm1 <- subset(tpm, rownames(tpm) %in% rownames(mean_tpm_keep)) #keep the top 75% most variable transcripts
dim(tpm1)[1]
tpm1 <- tpm1[which(row.names(tpm1) %in% coding$V1),] #keep only coding transcripts
tpm1[] <- lapply(tpm1, function(x) as.numeric(as.character(x)))

#center and scale the data, for further analysis
tpm2 <- apply(tpm1, 1, function(x) (x-mean(x))/sd(x))
tpm2 <- t(tpm2)
tpm2 <- as.data.frame(tpm2)
dim(tpm2)[1] #number of transcripts for downstream

###hierarchical clustering =====
hc <- hclust(dist(tpm2, method = "euclidean"), "complete")
save(hc, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/hierClust.rda")
dend <- as.dendrogram(hc)
save(dend, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/dendrogram.rda")

#cutree to make clusters
trans_cluster <- cutree(hc, k = 3) %>% enframe()
table(trans_cluster$value)

###summarize the time points by mean scaled tpm ====
tpm3 <- tpm2
tpm3 <- cbind(rownames(tpm2), tpm3)
rownames(tpm3) <- 1:nrow(tpm3)
colnames(tpm3)[1] <- "transcripts"

E7.5_mean <- rowMeans(tpm3[,c("E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
E8.5_mean <- rowMeans(tpm3[,c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
E9.5_mean <- rowMeans(tpm3[,c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5")])

#e7.5
e7.5_table <- data.frame(tpm3$transcripts, E7.5_mean)
e7.5_table$time <- c("e7.5")
rownames(e7.5_table) <- 1:nrow(e7.5_table)
colnames(e7.5_table) <- c("transcripts", "mean_cts_scaled", "time")

#e8.5
e8.5_table <- data.frame(tpm3$transcripts, E8.5_mean)
e8.5_table$time <- c("e8.5")
rownames(e8.5_table) <- 1:nrow(e8.5_table)
colnames(e8.5_table) <- c("transcripts", "mean_cts_scaled", "time")

#e9.5
e9.5_table <- data.frame(tpm3$transcripts, E9.5_mean)
e9.5_table$time <- c("e9.5")
rownames(e9.5_table) <- 1:nrow(e9.5_table)
colnames(e9.5_table) <- c("transcripts", "mean_cts_scaled", "time")

summary <- rbind(rbind(e7.5_table, e8.5_table), e9.5_table)
summary <- summary[order(summary$transcripts),]
row.names(summary) <- 1:nrow(summary)
summary <- inner_join(summary, trans_cluster, by = c("transcripts" = "name"))
summary <- inner_join(summary, t2g, by = c("transcripts" = "target_id"))
write.table(summary, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/transcriptGroups.txt", sep  = "\t", row.names = F, quote = F)

#group 1
geneGroup1 <- subset(summary, summary$value == "1")
length(unique(geneGroup1$transcripts))
length(unique(geneGroup1$ens_gene))
write.table(geneGroup1, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/e8.5Group.txt", sep  = " ", quote = F)

#group 2
geneGroup2 <- subset(summary, summary$value == "2")
length(unique(geneGroup2$transcripts))
length(unique(geneGroup2$ens_gene))
write.table(geneGroup2, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/e9.5Group.txt", sep  = " ", quote = F)

#group 3
geneGroup3 <- subset(summary, summary$value == "3")
length(unique(geneGroup2$transcripts))
length(unique(geneGroup3$ens_gene))
write.table(geneGroup3, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2B_hierarchicalClustering/e7.5Group.txt", sep  = " ", quote = F)