#GO term plot:
library(ggplot2)

#e7.5 BP
greatE7.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_BP_filtered.csv", header = T, sep = ",")
greatE7.5 <- greatE7.5[order(greatE7.5$Hyper.Qval),]
greatE7.5 <- greatE7.5[c(which(greatE7.5$Terms == "placenta development"), grep("hypoxia|trophoblast", greatE7.5$Terms)),]

pdf("Z:/hhvu/Project1_2/IFPA/e7.5specific_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE7.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#A5B557", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E7.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
dev.off()

#e7.5 singleKO
greatE7.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e7.5specific_singleKO_all.tsv", header = F, sep = "\t", quote = "")
greatE7.5 <- subset(greatE7.5, greatE7.5$V3 <= 0.05 & greatE7.5$V4 > 2)
greatE7.5 <- greatE7.5[, -8]
colnames(greatE7.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE7.5 <- greatE7.5[order(greatE7.5$Hyper.Qval),]

greatE7.5 <- greatE7.5[c(which(greatE7.5$Terms == "abnormal trophoblast layer morphology"), grep("abnormal mural trophectoderm morphology|abnormal wound healing", greatE7.5$Terms)),]

pdf("Z:/hhvu/Project1_2/IFPA/e7.5specific_importantSingleKOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE7.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#A5B557", width=0.4) + coord_flip() +
  labs(title = "Interesting MGI Single KO Mouse Phenotype Terms \n E7.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
dev.off()

#e8.5 BP
greatE8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e8.5specific_BP_all.tsv", header = F, sep = "\t", quote = "")
greatE8.5 <- subset(greatE8.5, greatE8.5$V3 <= 0.05 & greatE8.5$V4 > 2)
greatE8.5 <- greatE8.5[, -8]
colnames(greatE8.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE8.5 <- greatE8.5[order(greatE8.5$Hyper.Qval),]

greatE8.5 <- greatE8.5[grep("embryonic placenta|labyrinthine", greatE8.5$Terms),]

pdf("Z:/hhvu/Project1_2/IFPA/e8.5specific_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE8.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#354E71", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E8.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
dev.off()

#e8.5 singleKO
greatE8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e8.5specific_singleKO_all.tsv", header = F, sep = "\t", quote = "")
greatE8.5 <- subset(greatE8.5, greatE8.5$V3 <= 0.05 & greatE8.5$V4 > 2)
greatE8.5 <- greatE8.5[, -8]
colnames(greatE8.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE8.5 <- greatE8.5[order(greatE8.5$Hyper.Qval),]

greatE8.5 <- greatE8.5[c(which(greatE8.5$Terms == "embryonic lethality"), grep("failure of chorioallantoic fusion|abnormal fetal growth/weight/body size", greatE8.5$Terms)),]

pdf("Z:/hhvu/Project1_2/IFPA/e8.5specific_importantSingleKOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE8.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#354E71", width=0.4) + coord_flip() +
  labs(title = "Interesting MGI Single KO Mouse Phenotype Terms \n E8.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=18), axis.title.x = element_text(size = 18))
dev.off()

#e9.5 BP
greatE9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e9.5specific_BP_all.tsv", header = F, sep = "\t", quote = "")
greatE9.5 <- subset(greatE9.5, greatE9.5$V3 <= 0.05 & greatE9.5$V4 > 2)
greatE9.5 <- greatE9.5[, -8]
colnames(greatE9.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE9.5 <- greatE9.5[order(greatE9.5$Hyper.Qval),]

greatE9.5 <- greatE9.5[grep("^regulation of angiogenesis|labyrinthine layer blood vessel development|regulation of cellular response to insulin stimulus", greatE9.5$Terms),]

pdf("Z:/hhvu/Project1_2/IFPA/e9.5specific_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE9.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#841F27", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E9.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
dev.off()

#e9.5 singleKO
greatE9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e9.5specific_singleKO_all.tsv", header = F, sep = "\t", quote = "")
greatE9.5 <- subset(greatE9.5, greatE9.5$V3 <= 0.05 & greatE9.5$V4 > 2)
greatE9.5 <- greatE9.5[, -8]
colnames(greatE9.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE9.5 <- greatE9.5[order(greatE9.5$Hyper.Qval),]

greatE9.5 <- greatE9.5[c(which(greatE9.5$Terms == "abnormal placenta labyrinth morphology"), grep("abnormal placenta labyrinth size|decreased total body fat amount", greatE9.5$Terms)),]

pdf("Z:/hhvu/Project1_2/IFPA/e9.5specific_importantSingleKOTerms.pdf", width = 10, height = 5)
ggplot(data=greatE9.5, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#841F27", width=0.4) + coord_flip() +
  labs(title = "Interesting MGI Single KO Mouse Phenotype Terms \n E9.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=18), axis.title.x = element_text(size = 18))
dev.off()