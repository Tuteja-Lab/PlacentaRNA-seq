library("ggplot2")

#cluster profiler for subnetworks - plots:

#e7.5 BP
load("Z:/hhvu/Project1_2/RNA-seq/6_clusterProfiler/GENIE3/max1000/e7.5_2_BP_filtered.rda")
goTerms <- goTerms[c(which(goTerms$Description == "epithelial cell proliferation"), which(goTerms$Description == "epithelium migration"), which(goTerms$Description == "cell differentiation involved in embryonic placenta development")),]
goTerms$Description[3] <- "cell differentiation involved in \n embryonic placenta development"

pdf("Z:/hhvu/Project1_2/IFPA/e7.5_2_GENIE3_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=goTerms, aes(x=Description, y=-log10(as.numeric(qvalue)))) +
  geom_bar(stat="identity", fill="#A5B557", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E7.5 Subnetwork - GENIE3",
       y = "-log10(q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
dev.off()

#e8.5 BP
load("Z:/hhvu/Project1_2/RNA-seq/6_clusterProfiler/STRING/max1000/e8.5_1_BP_filtered.rda")
goTerms <- goTerms[c(which(goTerms$Description == "tube morphogenesis"), which(goTerms$Description == "placenta development"), 
                     which(goTerms$Description == "vasculature development")),]

pdf("Z:/hhvu/Project1_2/IFPA/e8.5_1_STRING_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=goTerms, aes(x=Description, y=-log10(as.numeric(qvalue)))) +
  geom_bar(stat="identity", fill="#354E71", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E8.5 Subnetwork - STRING",
       y = "-log10(q-value)",
       x = "") +
  theme(plot.title = element_text(size = 20, face = "bold"), legend.text=element_text(size=20),
        axis.text=element_text(size=20), axis.title.x = element_text(size = 20))
dev.off()

#e9.5 BP
load("Z:/hhvu/Project1_2/RNA-seq/6_clusterProfiler/GENIE3/max1000/e9.5_3_BP_filtered.rda")
goTerms <- goTerms[c(which(goTerms$Description == "blood vessel development"), which(goTerms$Description == "placenta development"), 
                     which(goTerms$Description == "endothelial cell proliferation")),]

pdf("Z:/hhvu/Project1_2/IFPA/e9.5_3_GENIE3_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=goTerms, aes(x=Description, y=-log10(as.numeric(qvalue)))) +
  geom_bar(stat="identity", fill="#841F27", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E9.5 Subnetwork - STRING",
       y = "-log10(q-value)",
       x = "") +
  theme(plot.title = element_text(size = 20, face = "bold"), legend.text=element_text(size=20),
        axis.text=element_text(size=20), axis.title.x = element_text(size = 20))
dev.off()

#GREAT for time-point specific genes
#e8.5 BP
greatE8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/4_specificGenes/GREAT/e8.5specific_BP_all.tsv", header = F, sep = "\t", quote = "")
greatE8.5 <- subset(greatE8.5, greatE8.5$V3 <= 0.05 & greatE8.5$V4 > 2)
greatE8.5 <- greatE8.5[, -8]
colnames(greatE8.5) <- c("Terms", "Hyper.Rank", "Hyper.Qval", "Fold", "Hyper.Obs.Gene.Hits", "Hyper.Total.Genes", "Hyper.Gene.Set.Cov")
greatE8.5 <- greatE8.5[order(greatE8.5$Hyper.Qval),]

greatE8.5 <- greatE8.5[grep("embryonic placenta|labyrinthine", greatE8.5$Terms),]

pdf("Z:/hhvu/Project1_2/IFPA/e8.5specific_importantGOTerms.pdf", width = 10, height = 5)
ggplot(data=goTerms, aes(x=Terms, y=-log10(as.numeric(Hyper.Qval)))) +
  geom_bar(stat="identity", fill="#354E71", width=0.4) + coord_flip() +
  labs(title = "Interesting GO Biological Process Terms \n E8.5-specific Genes",
       y = "-log10(Hypergeomatric q-value)",
       x = "") +
  theme(plot.title = element_text(size = 15, face = "bold"), legend.text=element_text(size=15),
        axis.text=element_text(size=15), axis.title.x = element_text(size = 15))
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
