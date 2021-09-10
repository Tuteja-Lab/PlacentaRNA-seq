#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("devtools")    # only if devtools not yet installed
#BiocManager::install("pachterlab/sleuth")

library("sleuth")
library("dplyr")
library("biomaRt")

set.seed(123)

###building sample description
dir <- "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/1_kallisto/official/"
sample_id <- dir(file.path(dir))
sample_id <- sample_id[!(sample_id %in% c("estCounts_allSamples.tsv", "TPM_allSamples.tsv"))] #remove the unrelated file if there is

kal_dirs <- file.path(dir, sample_id, "abundance.h5")

desc <- matrix(c(sample_id, "e9.5", "e8.5", "e9.5", "e9.5", "e9.5", "e9.5", "e9.5", "e7.5", "e7.5", "e7.5", "e7.5", "e7.5", "e7.5", "e8.5", "e8.5", "e8.5", "e8.5", "e8.5"), nrow = 18, ncol = 2, byrow = FALSE)
desc <- as.data.frame(desc)
desc <- dplyr::mutate(desc, path = kal_dirs)
colnames(desc) <- c("sample", "condition", "path")
desc <- desc[order(desc$condition), ]
desc <- subset(desc, !(desc$sample %in% c("S55", "S56"))) #remove sample S55, S56 outlier

e7.5vsE8.5 <- subset(desc, (desc$condition %in% c("e7.5", "e8.5")))
e8.5vsE9.5 <- subset(desc, (desc$condition %in% c("e9.5", "e8.5")))
e7.5vsE9.5 <- subset(desc, (desc$condition %in% c("e7.5", "e9.5")))

###building transcript name table. Since the object is presaved, just load it here
#still include the code to build t2g object
#mart <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL', dataset = "mmusculus_gene_ensembl", host = "http://sep2019.archive.ensembl.org") 
#t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "external_gene_name"), mart = mart)
#t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id_version, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

t2g <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/t2g.txt", header = T, sep = "\t")

#load protein-coding transcripts
coding <- read.table("/work/LAS/geetu-lab/hhvu/Mus_musculus_grcm38_coding_transcripts.txt", header = F)

###function to calculate sleuth genes + transcripts, likelihood ratio test

#this function will return a result table (named deg) with DE protein-coding transcripts that belong to DE genes
#and a table of all transcripts + necessary metrics (named trans), for plotting
deg <-  function(description, targetMap, con1, con2){
  #### 1.1. gene level
  so <- sleuth_prep(description, target_mapping = targetMap, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
  so <- sleuth_lrt(so, 'reduced', 'full')
  
  sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = so$pval_aggregate)
  sleuth_table <- dplyr::filter(sleuth_table, qval <= 0.05)
  
  #### 1.2. transcript level
  soTrans <- sleuth_prep(description, target_mapping = targetMap, extra_bootstrap_summary = TRUE)
  soTrans <- sleuth_fit(soTrans, ~condition, 'full')
  soTrans <- sleuth_fit(soTrans, ~1, 'reduced')
  soTrans <- sleuth_lrt(soTrans, 'reduced', 'full')
  
  sleuth_tableTrans <- sleuth_results(soTrans, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = soTrans$pval_aggregate)
  
  ## finding up/down regulated transcripts
  cond1 <- description$sample[description$condition == con1]
  cond2 <- description$sample[description$condition == con2]
  
  raw_cond1 <- subset(soTrans$obs_raw, soTrans$obs_raw$sample %in% cond1)
  raw_cond2 <- subset(soTrans$obs_raw, soTrans$obs_raw$sample %in% cond2)
  
  #calculate mean TPM of each transcript across all samples
  means_cond1 <- aggregate(tpm~target_id, data = raw_cond1, FUN=function(x) c(mean=mean(x)))
  colnames(means_cond1) <- c("target_id", "TPM_con1")
  
  means_cond2 <- aggregate(tpm~target_id, data = raw_cond2, FUN=function(x) c(mean=mean(x)))
  colnames(means_cond2) <- c("target_id", "TPM_con2")
  
  merged_means <- merge(means_cond1, means_cond2, by = c("target_id"))
  
  sleuth_tableTrans_TPM <- left_join(sleuth_tableTrans, merged_means)
  
  #calculate FC between 2 condition
  sleuth_tableTrans_TPM$TPM_con1 <- sleuth_tableTrans_TPM$TPM_con1 + 1
  sleuth_tableTrans_TPM$TPM_con2 <- sleuth_tableTrans_TPM$TPM_con2 + 1
  sleuth_tableTrans_TPM <- transform(sleuth_tableTrans_TPM, log2FC = log2(sleuth_tableTrans_TPM$TPM_con2 / sleuth_tableTrans_TPM$TPM_con1))
  
  #keep only DE transcripts of DE genes
  sleuth_tableTrans_TPM <- subset(sleuth_tableTrans_TPM,
                                  sleuth_tableTrans_TPM$ens_gene %in% sleuth_table$target_id
                                  & sleuth_tableTrans_TPM$qval <= 0.05
                                  & abs(sleuth_tableTrans_TPM$log2FC) >= log2(1.5))
  sleuth_tableTrans_TPM <- sleuth_tableTrans_TPM[order(sleuth_tableTrans_TPM$ext_gene),]
  
  #filter out non coding transcripts
  sleuth_tableTrans_TPM <- sleuth_tableTrans_TPM[sleuth_tableTrans_TPM$target_id %in% coding$V1,]
  sleuth_tableTrans <- sleuth_tableTrans[sleuth_tableTrans$target_id %in% coding$V1,]
  
  
  results <- list("deg" = sleuth_tableTrans_TPM, "trans" = sleuth_tableTrans)
  return(results)
}


### transcript level + gene level. Only get DE trans belonging to DE genes: #####
#Notes: There may be some warnings like this:
#"NA values were found during variance shrinkage estimation due to mean observation values outside of the range used for the LOESS fit"
#Since the number of NA transcripts is small (only 1 to 6 transcripts), and the NA values won't affect the adjustment of p-values
#we can ignore the warnings.

# e7.5 vs e8.5 test
e7.5_e8.5 <- deg(e7.5vsE8.5, t2g, "e7.5", "e8.5")

e7.5vsE8.5deg <- e7.5_e8.5$deg
length(unique(e7.5vsE8.5deg$target_id)) #check total number of DE transcripts
length(unique(e7.5vsE8.5deg$ens_gene)) #check total number of DE genes

e7.5_e7.5vsE8.5 <- e7.5vsE8.5deg[e7.5vsE8.5deg$log2FC <= -log2(1.5),]
length(unique(e7.5_e7.5vsE8.5$target_id)) #check total number of transcripts up-regulated at e7.5, compared to e8.5
length(unique(e7.5_e7.5vsE8.5$ens_gene)) #check total number of genes up-regulated at e7.5, compared to e8.5
write.table(e7.5_e7.5vsE8.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e7.5_e7.5vsE8.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")

e8.5_e7.5vsE8.5 <- e7.5vsE8.5deg[e7.5vsE8.5deg$log2FC >= log2(1.5),]
length(unique(e8.5_e7.5vsE8.5$target_id)) #check total number of transcripts up-regulated at e8.5, compared to e7.5
length(unique(e8.5_e7.5vsE8.5$ens_gene)) #check total number of genes up-regulated at e8.5, compared to e7.5
write.table(e8.5_e7.5vsE8.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e8.5_e7.5vsE8.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")


# e7.5 vs e9.5 test
e7.5_e9.5 <- deg(e7.5vsE9.5, t2g, "e7.5", "e9.5")

e7.5vse9.5deg <- e7.5_e9.5$deg
length(unique(e7.5vse9.5deg$target_id)) #check total number of DE transcripts
length(unique(e7.5vse9.5deg$ens_gene)) #check total number of DE genes

e7.5_e7.5vse9.5 <- e7.5vse9.5deg[e7.5vse9.5deg$log2FC <= -log2(1.5),]
length(unique(e7.5_e7.5vse9.5$target_id)) #check total number of transcripts up-regulated at e7.5, compared to e9.5
length(unique(e7.5_e7.5vse9.5$ens_gene)) #check total number of genes up-regulated at e7.5, compared to e9.5
write.table(e7.5_e7.5vse9.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e7.5_e7.5vse9.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")

e9.5_e7.5vse9.5 <- e7.5vse9.5deg[e7.5vse9.5deg$log2FC >= log2(1.5),]
length(unique(e9.5_e7.5vse9.5$target_id)) #check total number of transcripts up-regulated at e9.5, compared to e7.5
length(unique(e9.5_e7.5vse9.5$ens_gene)) #check total number of genes up-regulated at e9.5, compared to e7.5
write.table(e9.5_e7.5vse9.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e9.5_e7.5vse9.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")


# e8.5 vs e9.5 test
e8.5_e9.5 <- deg(e8.5vsE9.5, t2g, "e8.5", "e9.5")

e8.5vse9.5deg <- e8.5_e9.5$deg
length(unique(e8.5vse9.5deg$target_id)) #check total number of DE transcripts
length(unique(e8.5vse9.5deg$ens_gene)) #check total number of DE genes

e8.5_e8.5vse9.5 <- e8.5vse9.5deg[e8.5vse9.5deg$log2FC <= -log2(1.5),]
length(unique(e8.5_e8.5vse9.5$target_id)) #check total number of transcripts up-regulated at e8.5, compared to e9.5
length(unique(e8.5_e8.5vse9.5$ens_gene)) #check total number of genes up-regulated at e8.5, compared to e9.5
write.table(e8.5_e8.5vse9.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e8.5_e8.5vse9.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")

e9.5_e8.5vse9.5 <- e8.5vse9.5deg[e8.5vse9.5deg$log2FC >= log2(1.5),]
length(unique(e9.5_e8.5vse9.5$target_id)) #check total number of transcripts up-regulated at e9.5, compared to e8.5
length(unique(e9.5_e8.5vse9.5$ens_gene)) #check total number of genes up-regulated at e9.5, compared to e8.5
write.table(e9.5_e8.5vse9.5, "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e9.5_e8.5vse9.5_DEtransGenes.txt", row.names = F, quote = F, sep = "\t")


#to do: volcano plots (20210701)
#function deg2 is for plotting purpose, as for volcano plot I needed insignificant transcripts also
deg2 <-  function(description, targetMap, con1, con2){
  #### 1.2. transcript level
  so1trans <- sleuth_prep(description, target_mapping = targetMap, extra_bootstrap_summary = TRUE)
  so1trans <- sleuth_fit(so1trans, ~condition, 'full')
  so1trans <- sleuth_fit(so1trans, ~1, 'reduced')
  so1trans <- sleuth_lrt(so1trans, 'reduced', 'full')
  
  sleuth_table1trans <- sleuth_results(so1trans, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = so1trans$pval_aggregate)
  
  ## finding up/down regulated transcripts
  cond1 <- description$sample[description$condition == con1]
  cond2 <- description$sample[description$condition == con2]
  
  raw_cond1 <- subset(so1trans$obs_raw, so1trans$obs_raw$sample %in% cond1)
  raw_cond2 <- subset(so1trans$obs_raw, so1trans$obs_raw$sample %in% cond2)
  
  #calculate mean TPM of each transcript across all samples
  means_cond1 <- aggregate(tpm~target_id, data = raw_cond1, FUN=function(x) c(mean=mean(x)))
  colnames(means_cond1) <- c("target_id", "TPM_con1")
  
  means_cond2 <- aggregate(tpm~target_id, data = raw_cond2, FUN=function(x) c(mean=mean(x)))
  colnames(means_cond2) <- c("target_id", "TPM_con2")
  
  merged_means <- merge(means_cond1, means_cond2, by = c("target_id"))
  
  sleuth_tableTrans_TPM <- left_join(sleuth_table1trans, merged_means)
  
  #calculate FC between 2 condition
  sleuth_tableTrans_TPM$TPM_con1 <- sleuth_tableTrans_TPM$TPM_con1 + 1
  sleuth_tableTrans_TPM$TPM_con2 <- sleuth_tableTrans_TPM$TPM_con2 + 1
  sleuth_tableTrans_TPM <- transform(sleuth_tableTrans_TPM, log2FC = log2(sleuth_tableTrans_TPM$TPM_con2 / sleuth_tableTrans_TPM$TPM_con1))
  
  return(sleuth_tableTrans_TPM)
}

e7.5_e8.5 <- deg2(e7.5vsE8.5, t2g, "e7.5", "e8.5")
save(e7.5_e8.5, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e7.5_e8.5_allTrans-forPlots.rda")

e7.5_e9.5 <- deg2(e7.5vsE9.5, t2g, "e7.5", "e9.5")
save(e7.5_e9.5, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e7.5_e9.5_allTrans-forPlots.rda")

e8.5_e9.5 <- deg2(e8.5vsE9.5, t2g, "e8.5", "e9.5")
save(e8.5_e9.5, file = "/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/2A_sleuth/e8.5_e9.5_allTrans-forPlots.rda")