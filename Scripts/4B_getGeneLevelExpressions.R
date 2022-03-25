library("tximport")
library("GENIE3")

###building transcript names
t2g <- read.table("Z:/hhvu/Project1_2/RNA-seq/t2g.txt", header = T, sep = "\t")

###import count files
sample_id <- dir(file.path("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/"))
sample_id <- sample_id[!(sample_id %in% c("estCounts_allSamples.tsv", "TPM_allSamples.tsv", "e7.5geneLevelTPM.rda",     "e8.5geneLevelTPM.rda",     "e9.5geneLevelTPM.rda"))] #remove the unrelated file if there is

kal_dirs <- file.path("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/", sample_id, "abundance.tsv")

desc <- matrix(c(sample_id, "e9.5", "e8.5", "e9.5", "e9.5", "e9.5", "e9.5", "e9.5", "e7.5", "e7.5", "e7.5", "e7.5", "e7.5", "e7.5", "e8.5", "e8.5", "e8.5", "e8.5", "e8.5"), nrow = 18, ncol = 2, byrow = FALSE)
desc <- as.data.frame(desc)
desc <- dplyr::mutate(desc, path = kal_dirs)
colnames(desc) <- c("sample", "condition", "path")
desc <- desc[order(desc$condition), ]
desc <- subset(desc, !(desc$sample %in% c("S55", "S56"))) #remove sample S55, S56 outlier


tpm <- read.table("Z:/hhvu/Project1_2/RNA-seq/1_kallisto/official/TPM_allSamples.tsv", header = T)
tpm <- tpm[order(tpm$target_id),]
row.names(tpm) <- tpm$target_id
tpm <- tpm[,-c(1, 20, 21)]
tpm <- tpm[, -which(names(tpm) %in% c("E7.5_1", "E9.5_6"))]
tpm[] <- lapply(tpm, function(x) as.numeric(as.character(x)))

#e7.5
e7.5Files <- subset(desc, (desc$condition %in% c("e7.5")))
e7.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e7.5specific_trans.txt", header = F)

keep <- subset(tpm, row.names(tpm) %in% e7.5$V1)
keep$E7.5_mean <- rowMeans(keep[,c("E7.5_2", "E7.5_3", "E7.5_4", "E7.5_5", "E7.5_6")])
keep <- subset(keep, keep$E7.5_mean >= 5)

e7.5_2 <- read.table(e7.5Files$path[1], header = T)
e7.5_2 <- subset(e7.5_2, e7.5_2$target_id %in% row.names(keep))
write.table(e7.5_2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e7.5_2.tsv", row.names = F, quote = F, sep = '\t')

e7.5_3 <- read.table(e7.5Files$path[2], header = T)
e7.5_3 <- subset(e7.5_3, e7.5_3$target_id %in% row.names(keep))
write.table(e7.5_3, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e7.5_3.tsv", row.names = F, quote = F, sep = '\t')

e7.5_4 <- read.table(e7.5Files$path[3], header = T)
e7.5_4 <- subset(e7.5_4, e7.5_4$target_id %in% row.names(keep))
write.table(e7.5_4, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e7.5_4.tsv", row.names = F, quote = F, sep = '\t')

e7.5_5 <- read.table(e7.5Files$path[4], header = T)
e7.5_5 <- subset(e7.5_5, e7.5_5$target_id %in% row.names(keep))
write.table(e7.5_5, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e7.5_5.tsv", row.names = F, quote = F, sep = '\t')

e7.5_6 <- read.table(e7.5Files$path[5], header = T)
e7.5_6 <- subset(e7.5_6, e7.5_6$target_id %in% row.names(keep))
write.table(e7.5_6, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e7.5_6.tsv", row.names = F, quote = F, sep = '\t')


#e8.5
e8.5Files <- subset(desc, (desc$condition %in% c("e8.5")))
e8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e8.5specific_trans.txt", header = F)

keep <- subset(tpm, row.names(tpm) %in% e8.5$V1)
keep$E8.5_mean <- rowMeans(keep[,c("E8.5_1", "E8.5_2", "E8.5_3", "E8.5_4", "E8.5_5", "E8.5_6")])
keep <- subset(keep, keep$E8.5_mean >= 5)

e8.5_1 <- read.table(e8.5Files$path[1], header = T)
e8.5_1 <- subset(e8.5_1, e8.5_1$target_id %in% row.names(keep))
write.table(e8.5_1, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_1.tsv", row.names = F, quote = F, sep = '\t')

e8.5_2 <- read.table(e8.5Files$path[2], header = T)
e8.5_2 <- subset(e8.5_2, e8.5_2$target_id %in% row.names(keep))
write.table(e8.5_2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_2.tsv", row.names = F, quote = F, sep = '\t')

e8.5_3 <- read.table(e8.5Files$path[3], header = T)
e8.5_3 <- subset(e8.5_3, e8.5_3$target_id %in% row.names(keep))
write.table(e8.5_3, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_3.tsv", row.names = F, quote = F, sep = '\t')

e8.5_4 <- read.table(e8.5Files$path[4], header = T)
e8.5_4 <- subset(e8.5_4, e8.5_4$target_id %in% row.names(keep))
write.table(e8.5_4, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_4.tsv", row.names = F, quote = F, sep = '\t')

e8.5_5 <- read.table(e8.5Files$path[5], header = T)
e8.5_5 <- subset(e8.5_5, e8.5_5$target_id %in% row.names(keep))
write.table(e8.5_5, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_5.tsv", row.names = F, quote = F, sep = '\t')

e8.5_6 <- read.table(e8.5Files$path[6], header = T)
e8.5_6 <- subset(e8.5_6, e8.5_6$target_id %in% row.names(keep))
write.table(e8.5_6, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e8.5_6.tsv", row.names = F, quote = F, sep = '\t')


#e9.5
e9.5Files <- subset(desc, (desc$condition %in% c("e9.5")))
e9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/e9.5specific_trans.txt", header = F)

keep <- subset(tpm, row.names(tpm) %in% e9.5$V1)
keep$e9.5_mean <- rowMeans(keep[,c("E9.5_1", "E9.5_2", "E9.5_3", "E9.5_4", "E9.5_5")])
keep <- subset(keep, keep$e9.5_mean >= 5)

e9.5_1 <- read.table(e9.5Files$path[1], header = T)
e9.5_1 <- subset(e9.5_1, e9.5_1$target_id %in% row.names(keep))
write.table(e9.5_1, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e9.5_1.tsv", row.names = F, quote = F, sep = '\t')

e9.5_2 <- read.table(e9.5Files$path[2], header = T)
e9.5_2 <- subset(e9.5_2, e9.5_2$target_id %in% row.names(keep))
write.table(e9.5_2, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e9.5_2.tsv", row.names = F, quote = F, sep = '\t')

e9.5_3 <- read.table(e9.5Files$path[3], header = T)
e9.5_3 <- subset(e9.5_3, e9.5_3$target_id %in% row.names(keep))
write.table(e9.5_3, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e9.5_3.tsv", row.names = F, quote = F, sep = '\t')

e9.5_4 <- read.table(e9.5Files$path[4], header = T)
e9.5_4 <- subset(e9.5_4, e9.5_4$target_id %in% row.names(keep))
write.table(e9.5_4, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e9.5_4.tsv", row.names = F, quote = F, sep = '\t')

e9.5_5 <- read.table(e9.5Files$path[5], header = T)
e9.5_5 <- subset(e9.5_5, e9.5_5$target_id %in% row.names(keep))
write.table(e9.5_5, "Z:/hhvu/Project1_2/RNA-seq/5B_GENIE3/kallisto_GENIE3/e9.5_5.tsv", row.names = F, quote = F, sep = '\t')
