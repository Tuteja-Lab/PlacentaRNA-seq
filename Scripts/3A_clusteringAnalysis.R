library("dplyr")
library("tidyverse")
library("biomaRt")
library("ggplot2")
library("dendextend")
library("matrixStats")
library("RclusTool")
library("kohonen")

set.seed(123)

###building transcript names =====
t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/t2g.txt", header = T, sep = "\t")
t2g <- t2g[order(t2g$target_id),]
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
mean_tpm$var <- matrixStats::rowVars(as.matrix(mean_tpm)) #calculate the variance
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

## som
som <- som(as.matrix(tpm2), grid = somgrid(3, 1, "rectangular"))
trans_cluster_som <- som$unit.classif %>% enframe()
trans_cluster_som$name <- row.names(tpm2)

## k means
km <- kmeans(tpm2, centers = 3)
trans_cluster_kmeans <- km$cluster %>% enframe()


## spectral clustering - this may take a long time, should run on HPC
sim <- computeGaussianSimilarity(tpm2, 1)
res <- spectralClustering(sim, K=3)
trans_cluster_sc <- data_frame(row.names(tpm2), res$label)
colnames(trans_cluster_sc) <- c("name", "value")

write.table(trans_cluster_sc, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/3_clusterAnalysis/trans_cluster_sc.txt", quote = F, row.names = F)


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

summary <- rbind(e7.5_table, e8.5_table)
summary <- rbind(summary, e9.5_table)
summary <- summary[order(summary$transcripts),]
row.names(summary) <- 1:nrow(summary)

#IMPORTANT: change trans_cluster_* according to the methods
summary <- inner_join(summary, trans_cluster_sc, by = c("transcripts" = "name"))
#summary <- inner_join(summary, trans_cluster_som, by = c("transcripts" = "name"))
#summary <- inner_join(summary, trans_cluster_kmeans, by = c("transcripts" = "name"))

summary <- inner_join(summary, t2g, by = c("transcripts" = "target_id"))

#with marker genes
mmp9 <- subset(summary, summary$ext_gene == "Mmp9") #e7.5
adm <- subset(summary, summary$ext_gene == "Adm") #e7.5
twist1 <- subset(summary, summary$ext_gene == "Twist1") #e8.5
tbx4 <- subset(summary, summary$ext_gene == "Tbx4") #e9.5
grb2 <- subset(summary, summary$ext_gene == "Grb2") #e8.5
itga4 <- subset(summary, summary$ext_gene == "Itga4")
gcm1 <- subset(summary, summary$ext_gene == "Gcm1") #e9.5

pdf("dir/sc_groups_markers_0.75quantile.pdf", width=12, height=8)
ggplot(aes(time, mean_cts_scaled), data = summary) +
  geom_line(aes(group = transcripts), alpha = 0.5, colour = "grey77") + #grey
  geom_line(stat = "summary", fun = median, size = 1.2, aes(group = 1, color = "Group median")) +
  labs(title = "Spectral clustering of transcripts",
       x = "Time point",
       y = "Scaled mean transcript counts", color = "Legend", linetype = "Legend") +
  theme(plot.title = element_text(size = 20, face = "bold"), legend.text=element_text(size=20)) +
  geom_line(data = mmp9, size = 1.2, aes(group = transcripts, color = "Mmp9", linetype = "Mmp9"), alpha = 1) + #e7.5
  geom_line(data = adm, size = 1.2, aes(group = transcripts, color = "Adm", linetype = "Adm" ), alpha = 1) + #e7.5
  geom_line(data = twist1, size = 1.2, aes(group = transcripts, color = "Twist1", linetype = "Twist1"), alpha = 1) + #e8.5
  geom_line(data = tbx4, size = 1.2, aes(group = transcripts, color = "Tbx4", linetype = "Tbx4"), alpha = 1) + #e9.5
  geom_line(data = itga4, size = 1.2, aes(group = transcripts, color = "Itga4", linetype = "Itga4"), alpha = 1) + #e8.5
  geom_line(data = gcm1, size = 1.2, aes(group = transcripts, color = "Gcm1", linetype = "Gcm1"), alpha = 1) + #e9.5
  scale_color_manual(name = "Legend", values = c("Mmp9" = "darkolivegreen4", "Adm" = "yellow3",
                                                 "Twist1" = "dodgerblue4", "Tbx4" = "saddlebrown",
                                                 "Itga4" = "deepskyblue4", "Gcm1" = "salmon3",
                                                 "Group median" = "grey22")) + 
  scale_linetype_manual(name = "Legend", values = c("Mmp9" = "solid", "Adm" = "solid",
                                                    "Twist1" = "twodash", "Tbx4" = "solid",
                                                    "Itga4" = "solid", "Gcm1" = "solid",
                                                    "Group median" = "solid")) + 
  guides(linetype=F,
         colour=guide_legend(keywidth = 3, keyheight = 1)) +
  theme(text = element_text(size=20), axis.text.x = element_text(angle=0, hjust=0.5)) +
  facet_grid(cols = vars(value))
dev.off()



