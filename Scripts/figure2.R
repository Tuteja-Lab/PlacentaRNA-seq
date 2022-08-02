library(clusterProfiler)
library("org.Mm.eg.db")
library(ggplot2)

go2 <- function(subnetworkGenes) {
  GO <- enrichGO(gene = subnetworkGenes, OrgDb=org.Mm.eg.db, ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, maxGSSize=1000, readable = T, keyType = "ENSEMBL")
  return(GO)
}

terms <- c("inflammatory response", "female pregnancy", "morphogenesis of a branching structure",
           "lipid biosynthetic process", "endothelial cell proliferation", "cholesterol metabolic process",
           "response to insulin", "placenta development", "vasculature development", 
           "positive regulation of cell migration", "epithelium migration")
terms2 <- c("inflammatory response", "female pregnancy", "morphogenesis of\na branching structure",
            "lipid biosynthetic process", "endothelial cell proliferation", "cholesterol metabolic process",
            "response to insulin", "placenta development", "vasculature development", 
            "positive regulation of cell migration", "epithelium migration")


#STRING
final <- data.frame(Description=as.character(),
                    GeneRatio=as.numeric(),
                    BgRatio=as.numeric(),
                    qvalue=as.numeric(),
                    Count=as.numeric(),
                    Fold=as.numeric(),
                    Rank=as.numeric(),
                    Rank2=as.numeric(),
                    Group=as.character())
type <- "STRING"
i <- 1
for (time in c("e7.5", "e8.5")) {
  geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/", time, "/largestComponent/", time, "_1/", time, "_1_TSS.bed"), header = F, sep = "\t")
  goTerms2 <- go2(geneList$V4)
  
  res <- as.data.frame(goTerms2@result)
  res$Rank <- seq(1, nrow(res))
  res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)
  
  res <- res[res$Description %in% terms,]
  res <- res[, c("Description", "GeneRatio", "BgRatio", "qvalue", "Count", "Fold", "Rank")]
  if (length(setdiff(terms, res$Description)) > 0) {
    ap <- data.frame(Description=setdiff(terms, res$Description),
                     GeneRatio=rep(0,length(setdiff(terms, res$Description))),
                     BgRatio=rep(0,length(setdiff(terms, res$Description))),
                     qvalue = rep(1, length(length(setdiff(terms, res$Description)))),
                     Count = rep(0, length(setdiff(terms, res$Description))),
                     Fold= rep(0, length(setdiff(terms, res$Description))),
                     Rank = rep(10000, length(setdiff(terms, res$Description))))
    res <- rbind(res, ap)
  }
  res <- res[order(res$qvalue),]
  res$Count <- as.numeric(res$Count)
  res$Rank2 <- ifelse(res$Rank <= 25, "Rank <= 25", 
                      ifelse(res$Rank > 25 & res$Rank <= 50, "25 < Rank <= 50",
                             ifelse(res$Rank > 50 & res$Rank <= 75, "50 < Rank <= 75",
                                    ifelse(res$Rank > 75 & res$Rank <= 100, "75 < Rank <= 100",
                                           ifelse(res$Rank > 100, "Rank > 100",
                                                  "black")))))
  res$Rank2 <- ifelse(res$qvalue > 0.05 | res$Fold < 2 | res$Count < 5, "Insignificant", res$Rank2) #
  res$Group <- paste0(time, "_", i, "_", type)
  res$Description <- ifelse(res$Description == "morphogenesis of a branching structure", "morphogenesis of\na branching structure", res$Description)
  res$Description <- factor(res$Description, levels = rev(terms2))
  final <- rbind(final, res)
}
time <- "e9.5"
for (i in 1:4) {
  geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5A_STRING/", time, "/largestComponent/", time, "_", i, "/", time, "_", i, "_TSS.bed"), header = F, sep = "\t")
  goTerms2 <- go2(geneList$V4)
  
  res <- as.data.frame(goTerms2@result)
  res$Rank <- seq(1, nrow(res))
  res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)
  #write.table(res[res$qvalue <= 0.05 & res$Fold >= 2 & res$Count >= 5,], 
  #            paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/", time, "_", i, "_", type, ".txt"), row.names = F, sep = "\t", quote = F)
  
  res <- res[res$Description %in% terms,]
  res <- res[, c("Description", "GeneRatio", "BgRatio", "qvalue", "Count", "Fold", "Rank")]
  
  if (length(setdiff(terms, res$Description)) > 0) {
    ap <- data.frame(Description=setdiff(terms, res$Description),
                     GeneRatio=rep(0,length(setdiff(terms, res$Description))),
                     BgRatio=rep(0,length(setdiff(terms, res$Description))),
                     qvalue = rep(1, length(length(setdiff(terms, res$Description)))),
                     Count = rep(0, length(setdiff(terms, res$Description))),
                     Fold = rep(0, length(setdiff(terms, res$Description))),
                     Rank = rep(10000, length(setdiff(terms, res$Description))))
    res <- rbind(res, ap)
  }
  res <- res[order(-res$qvalue),]
  res$Rank2 <- ifelse(res$Rank <= 25, "Rank <= 25", 
                      ifelse(res$Rank > 25 & res$Rank <= 50, "25 < Rank <= 50",
                             ifelse(res$Rank > 50 & res$Rank <= 75, "50 < Rank <= 75",
                                    ifelse(res$Rank > 75 & res$Rank <= 100, "75 < Rank <= 100",
                                           ifelse(res$Rank > 100, "Rank > 100",
                                                  "black")))))
  res$Rank2 <- ifelse(res$qvalue > 0.05 | res$Fold < 2 | res$Count < 5, "Insignificant", res$Rank2)
  res$Group <- paste0(time, "_", i, "_", type)
  res$Description <- ifelse(res$Description == "morphogenesis of a branching structure", "morphogenesis of\na branching structure", res$Description)
  res$Description <- factor(res$Description, levels = rev(terms2))
  final <- rbind(final, res)
}

#GENIE3
type <- "GENIE3"
i <- 2
for (time in c("e7.5", "e8.5")) {
  geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5B_GENIE3/", time, "/", time, "_", i, "/", time, "_", i, "_TSS.bed"), header = F, sep = "\t")
  goTerms2 <- go2(geneList$V4)
  
  res <- as.data.frame(goTerms2@result)
  res$Rank <- seq(1, nrow(res))
  res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)
  #write.table(res[res$qvalue <= 0.05 & res$Fold >= 2 & res$Count >= 5,], 
  #            paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/", time, "_", i, "_", type, ".txt"), row.names = F, sep = "\t", quote = F)
  
  res <- res[res$Description %in% terms,]
  res <- res[, c("Description", "GeneRatio", "BgRatio", "qvalue", "Count", "Fold", "Rank")]
  
  if (length(setdiff(terms, res$Description)) > 0) {
    ap <- data.frame(Description=setdiff(terms, res$Description),
                     GeneRatio=rep(0,length(setdiff(terms, res$Description))),
                     BgRatio=rep(0,length(setdiff(terms, res$Description))),
                     qvalue = rep(1, length(length(setdiff(terms, res$Description)))),
                     Count = rep(0, length(setdiff(terms, res$Description))),
                     Fold= rep(0, length(setdiff(terms, res$Description))),
                     Rank = rep(10000, length(setdiff(terms, res$Description))))
    res <- rbind(res, ap)
  }
  res <- res[order(-res$qvalue),]
  res$Rank2 <- ifelse(res$Rank <= 25, "Rank <= 25", 
                      ifelse(res$Rank > 25 & res$Rank <= 50, "25 < Rank <= 50",
                             ifelse(res$Rank > 50 & res$Rank <= 75, "50 < Rank <= 75",
                                    ifelse(res$Rank > 75 & res$Rank <= 100, "75 < Rank <= 100",
                                           ifelse(res$Rank > 100, "Rank > 100",
                                                  "black")))))
  res$Rank2 <- ifelse(res$qvalue > 0.05 | res$Fold < 2 | res$Count < 5, "Insignificant", res$Rank2)
  res$Group <- paste0(time, "_", i, "_", type)
  res$Description <- ifelse(res$Description == "morphogenesis of a branching structure", "morphogenesis of\na branching structure", res$Description)
  res$Description <- factor(res$Description, levels = rev(terms2))
  final <- rbind(final, res)
}
time <- "e9.5"
for (i in 1:3) {
  geneList <- read.table(paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/5B_GENIE3/", time, "/", time, "_", i, "/", time, "_", i, "_TSS.bed"), header = F, sep = "\t")
  goTerms2 <- go2(geneList$V4)
  
  res <- as.data.frame(goTerms2@result)
  res$Rank <- seq(1, nrow(res))
  res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)
  #write.table(res[res$qvalue <= 0.05 & res$Fold >= 2 & res$Count >= 5,], 
  #            paste0("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/", time, "_", i, "_", type, ".txt"), row.names = F, sep = "\t", quote = F)
  
  res <- res[res$Description %in% terms,]
  res <- res[, c("Description", "GeneRatio", "BgRatio", "qvalue", "Count", "Fold", "Rank")]
  
  if (length(setdiff(terms, res$Description)) > 0) {
    ap <- data.frame(Description=setdiff(terms, res$Description),
                     GeneRatio=rep(0,length(setdiff(terms, res$Description))),
                     BgRatio=rep(0,length(setdiff(terms, res$Description))),
                     qvalue = rep(1, length(length(setdiff(terms, res$Description)))),
                     Count = rep(0, length(setdiff(terms, res$Description))),
                     Fold= rep(0, length(setdiff(terms, res$Description))),
                     Rank = rep(10000, length(setdiff(terms, res$Description))))
    res <- rbind(res, ap)
  }
  res <- res[order(res$qvalue),]
  res$Rank2 <- ifelse(res$Rank <= 25, "Rank <= 25", 
                      ifelse(res$Rank > 25 & res$Rank <= 50, "25 < Rank <= 50",
                             ifelse(res$Rank > 50 & res$Rank <= 75, "50 < Rank <= 75",
                                    ifelse(res$Rank > 75 & res$Rank <= 100, "75 < Rank <= 100",
                                           ifelse(res$Rank > 100, "Rank > 100",
                                                  "black")))))
  res$Rank2 <- ifelse(res$qvalue > 0.05 | res$Fold < 2 | res$Count < 5, "Insignificant", res$Rank2)
  res$Group <- paste0(time, "_", i, "_", type)
  res$Description <- ifelse(res$Description == "morphogenesis of a branching structure", "morphogenesis of\na branching structure", res$Description)
  res$Description <- factor(res$Description, levels = rev(terms2))
  final <- rbind(final, res)
}

type<-c()
for (i in 1:nrow(final)) {
  j<-ifelse(length(grep("STRING", final$Group[i])) > 0, "STRING", "GENIE3")
  type <- c(type, j)
}
final$type <- type

col <- c("Rank <= 25" = "#7a0177", "25 < Rank <= 50" = "#c51b8a",
         "50 < Rank <= 75" = "#f768a1", "75 < Rank <= 100" = "#4dac26",
         "Rank > 100" = "#b8e186", "Insignificant" = "grey")
final$qvalue <- ifelse(final$qvalue < 0, 0, final$qvalue)
save(final, file="/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/6_clusterProfiler/plots/final.rda")


pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/fig2.pdf", width=56, height=25)
p <- ggplot(data=final, aes(x=log2(Fold), y=Description,
                            color=Rank2, size=-log10(as.numeric(qvalue)))) + 
  geom_point(show.legend = FALSE) + ylab("") + xlab("-Log10(q-value)") + scale_color_manual(values = col) +
  scale_size_continuous(range = c(2, 30)) +
  theme(legend.text = element_blank(), legend.title = element_blank(),
        axis.text.y = element_text(size = 55), strip.text = element_text(size = 40),
        axis.text.x = element_text(size = 30)) +
  facet_wrap(~Group, nrow=2)
p
dev.off()

pdf("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/fig2-legend.pdf", width=60, height=20)
p <- ggplot(data=final, aes(x=log2(Fold), y=Description,
                            color=Rank2, size=-log10(qvalue))) + 
  geom_point(show.legend = TRUE) + ylab("") + xlab("Log2(Fold change)") + scale_color_manual(values = col) +
  scale_size_continuous(range = c(2, 30)) +
  theme(legend.text = element_text(size = 30), legend.title = element_text(size = 30),
        axis.text.y = element_text(size = 20), strip.text = element_text(size = 30),
        axis.text.x = element_text(size = 20)) +
  facet_wrap(~Group, nrow=2)
p
dev.off()