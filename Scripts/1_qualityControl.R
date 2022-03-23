library("dplyr")
library("tidyverse")
library("biomaRt")
library("ggplot2")
library("dendextend")
library("matrixStats")
library("stats")
library("pca3d")
library("scatterplot3d")
library("pheatmap")

###filter based on raw est. counts - before filtering outliers =====
coding <- read.table("Z:/hhvu/Mus_musculus_grcm38_coding_transcripts.txt", header = F)

#load file
est_counts <- read.csv("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/estCounts_allSamples.tsv", header = TRUE, sep = "\t", stringsAsFactors=FALSE)

#reformat count file
est_counts <- est_counts[order(est_counts$target_id),]
rownames(est_counts) <- est_counts$target_id #make transcript names to be row names
est_counts <- est_counts[, -c(1, 20, 21)] #exclude non-count columns

#make mean count table
E7.5_mean <- rowMeans(est_counts[, c("E7.5_1", "E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
E8.5_mean <- rowMeans(est_counts[, c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
E9.5_mean <- rowMeans(est_counts[, c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5", "E9.5_6")])
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
keep <- apply(mean_cts, 1, basic_filter) #1 means apply the function to each row

###load TPM file in to do clustering =====
tpm <- read.csv("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/TPM_allSamples.tsv", header = TRUE, sep = "\t", stringsAsFactors=FALSE)

#reformat count file
tpm <- tpm[order(tpm$target_id),]
row.names(tpm) <- tpm$target_id
tpm <- tpm[,-c(1, 20, 21)]
tpm[] <- lapply(tpm, function(x) as.numeric(as.character(x)))

#filter out the low count transcripts according to the est_counts filtering above
tpm <- tpm[keep, ]

#calculate the variance between 3 time points
#here I take the average of each time point to make sure the variance is only calculated based on time points
#if not doing this, the variance may include biological variances between replicates also
E7.5_mean <- rowMeans(tpm[, c("E7.5_1", "E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
E8.5_mean <- rowMeans(tpm[, c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
E9.5_mean <- rowMeans(tpm[, c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5", "E9.5_6")])
mean_tpm <- data.frame(matrix(ncol = 3, nrow = nrow(tpm)))
row.names(mean_tpm) <- rownames(tpm)
colnames(mean_tpm) <- c("E7.5_mean", "E8.5_mean", "E9.5_mean")
mean_tpm$E7.5_mean <- E7.5_mean
mean_tpm$E8.5_mean <- E8.5_mean
mean_tpm$E9.5_mean <- E9.5_mean
mean_tpm$var <- rowVars(as.matrix(mean_tpm)) #calculate the variance
mean_tpm_keep <- subset(mean_tpm, mean_tpm$var > quantile(mean_tpm$var, 0.50)) #keep the top 50% most variable transcripts


#keep the top 50% most variable transcripts
tpm1 <- subset(tpm, rownames(tpm) %in% rownames(mean_tpm_keep))
tpm1[] <- lapply(tpm1, function(x) as.numeric(as.character(x)))

#center and scale the data, for further analysis
tpm2 <- apply(tpm1, 1, function(x) (x-mean(x))/sd(x))
tpm2 <- t(tpm2)
tpm2 <- as.data.frame(tpm2)


###do clustering by PCA =====
pca <- prcomp(t(tpm2))
summary(pca)
gr <- as.factor(c("e7.5", "e7.5", "e7.5", "e7.5", "e7.5", "e7.5",
                  "e8.5", "e8.5", "e8.5", "e8.5", "e8.5", "e8.5",
                  "e9.5", "e9.5", "e9.5", "e9.5", "e9.5", "e9.5")) #dumb way to group the samples, sorry...
mycols <- c("black", "red", "green")
names(mycols) <- c("e7.5", "e8.5", "e9.5")

#different plots
#3D plot
pdf("dir/PCA_allSamples.pdf", width=9, height=9)
scatterPlot <- with(pca, {
  s3d <- scatterplot3d(pca$x[,1:3], pch = 16, cex.symbols=2, cex.axis=1.5, cex.lab=2,
                       color = c("#A5B557", "#A5B557", "#A5B557",
                                 "#A5B557", "#A5B557", "#A5B557",
                                 "#354E71", "#354E71", "#354E71",
                                 "#354E71", "#354E71", "#354E71",
                                 "#841F27", "#841F27", "#841F27",
                                 "#841F27", "#841F27", "#841F27"),
                       angle = -40, box = F)
  s3d.coords <- s3d$xyz.convert(pca$x[,1], pca$x[,2], pca$x[,3])
  text(s3d.coords$x, s3d.coords$y, labels = ifelse(colnames(tpm1) == c("E7.5_1", "E9.5_6"), as.character(colnames(tpm1)), ''), pos=1, cex=1.5)
  }
  )
legend("topleft", pch = 16, legend = c("e7.5", "e8.5", "e9.5"), col = c("#A5B557", "#354E71", "#841F27"), pt.cex = 2, cex=1.5)
dev.off()

#2D PCA by ggplot
df_out <- as.data.frame(pca$x)
df_out$Group <- gr
#the following plot with PC2 and PC3 showed clearly E9.5_6 is an outlier
ggplot(df_out, aes(x=PC2, y=PC3, color=Group, label=ifelse(colnames(tpm1) == c("E7.5_1", "E9.5_6"), as.character(colnames(tpm1)), ''))) + 
  geom_text(hjust=0.3,size=5, show.legend = F) + geom_point(size=3) +
  labs(title="Principle component analysis", x="Principle component 2", y="Principle component 2") +
  theme(text = element_text(size=20), axis.text.x = element_text(angle=0, hjust=1))


###hierarchical clustering of samples - hclust ====
hc <- hclust(dist(t(tpm2)), "ave")
dend <- as.dendrogram(hc)
labels(dend)
labels_colors(dend) <- c("#A5B557", "#841F27", "#A5B557", "#A5B557", "#A5B557", "#A5B557", 
                         "#A5B557", "#354E71", "#354E71", "#354E71", "#354E71", "#354E71", 
                         "#354E71", "#841F27", "#841F27", "#841F27", "#841F27", "#841F27")

pdf("dir/hclust_allSamples.pdf", width=9, height=9)
dend %>% hang.dendrogram(hang = -1) %>% plot()
legend("topright", pch = 16, legend = c("e7.5", "e8.5", "e9.5"), col = c("#A5B557", "#354E71", "#841F27"), cex=1.5)
dev.off()