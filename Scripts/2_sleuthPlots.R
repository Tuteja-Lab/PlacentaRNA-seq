library(ggplot2)
library(ggrepel)

#Sleuth plots
volcanoPlot <- function(table, tp1, tp2, color1, color2) {
  names(table)[names(table) == 'qval'] <- 'padj'
  names(table)[names(table) == 'log2FC'] <- 'log2FoldChange'
  row.names(table) <- table$target_id
  table$color <- ifelse(table$sign == tp1, color1,
                        ifelse(table$sign == tp2, color2, "grey77"))
  
  return(table)
}

# volcano plots e7.5 vs e9.5
e7.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e7.5_e7.5vsE9.5_DEtransGenes.txt", sep = "\t", header = T)
e9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e9.5_e7.5vsE9.5_DEtransGenes.txt", sep = "\t", header = T)

load("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e7.5_e9.5_allTrans-forPlots.rda")
e7.5_e9.5_allTrans <- e7.5_e9.5
e7.5_e9.5_allTrans$sign <- ifelse(e7.5_e9.5_allTrans$target_id %in% e7.5$target_id, "e7.5",
                                  ifelse(e7.5_e9.5_allTrans$target_id %in% e9.5$target_id, "e9.5", 0))

e7.5_e9.5_allTrans <- volcanoPlot(e7.5_e9.5_allTrans, "e7.5", "e9.5", "#A5B557", "#841F27")

p1 <- ggplot(e7.5_e9.5_allTrans) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=color)) +
  ggtitle("E7.5 vs e9.5") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(adjusted p-value)") +
  scale_colour_identity("Legend", labels = c("E7.5 transcripts", "Not significant transcripts", "E9.5 transcripts"),
                        breaks = c("#A5B557", "grey77", "#841F27"), guide = "legend") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "black", size = 1) +
  theme(#legend.position = "right", 
    plot.title = element_text(size = 25, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))
pdf("Z:/hhvu/Project1_2/IFPA/e7.5_e9.5_volcano_RNAseq.pdf", width=15, height=9)
p1
dev.off()


#volcano plot e7.5 vs e8.5
e7.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e7.5_e7.5vse8.5_DEtransGenes.txt", sep = "\t", header = T)
e8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e8.5_e7.5vse8.5_DEtransGenes.txt", sep = "\t", header = T)

load("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e7.5_e8.5_allTrans-forPlots.rda")
e7.5_e8.5_allTrans <- e7.5_e8.5
e7.5_e8.5_allTrans$sign <- ifelse(e7.5_e8.5_allTrans$target_id %in% e7.5$target_id, "e7.5",
                                  ifelse(e7.5_e8.5_allTrans$target_id %in% e8.5$target_id, "e8.5", 0))

e7.5_e8.5_allTrans <- volcanoPlot(e7.5_e8.5_allTrans, "e7.5", "e8.5", "#A5B557", "#354E71")

p2 <- ggplot(e7.5_e8.5_allTrans) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=color)) +
  ggtitle("E7.5 vs e8.5") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(adjusted p-value)") +
  scale_colour_identity("Legend", labels = c("E7.5 transcripts", "Not significant transcripts", "E8.5 transcripts"),
                        breaks = c("#A5B557", "grey77", "#354E71"), guide = "legend") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "black", size = 1) +
  theme(#legend.position = "right", 
    plot.title = element_text(size = 25, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))

pdf("Z:/hhvu/Project1_2/IFPA/e7.5_e8.5_volcano_RNAseq.pdf", width=15, height=9)
p2
dev.off()

#volcano plot e8.5 vs e9.5
e8.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e8.5_e8.5vse9.5_DEtransGenes.txt", sep = "\t", header = T)
e9.5 <- read.table("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e9.5_e8.5vse9.5_DEtransGenes.txt", sep = "\t", header = T)

load("Z:/hhvu/Project1_2/RNA-seq/2A_sleuth/e8.5_e9.5_allTrans-forPlots.rda")
e9.5_e8.5_allTrans <- e8.5_e9.5
e9.5_e8.5_allTrans$sign <- ifelse(e9.5_e8.5_allTrans$target_id %in% e9.5$target_id, "e9.5",
                                  ifelse(e9.5_e8.5_allTrans$target_id %in% e8.5$target_id, "e8.5", 0))

e9.5_e8.5_allTrans <- volcanoPlot(e9.5_e8.5_allTrans, "e9.5", "e8.5", "#841F27", "#354E71")

p3 <- ggplot(e9.5_e8.5_allTrans) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=color)) +
  ggtitle("E8.5 vs e9.5") +
  xlab("log2(Fold Change)") + 
  ylab("-log10(adjusted p-value)") +
  scale_colour_identity("Legend", labels = c("E9.5 transcripts", "Not significant transcripts", "E8.5 transcripts"),
                        breaks = c("#841F27", "grey77", "#354E71"), guide = "legend") +
  #scale_y_continuous(limits = c(0,50)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = -log2(1.5), linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = log2(1.5), linetype = "dashed", color = "black", size = 1) +
  theme(#legend.position = "right", 
    plot.title = element_text(size = 25, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    legend.title = element_text(size=20),
    legend.text = element_text(size=20))

pdf("Z:/hhvu/Project1_2/IFPA/e8.5_e9.5_volcano_RNAseq.pdf", width=15, height=9)
p3
dev.off()

pdf("Z:/hhvu/Project1_2/IFPA/testGrid_volcano_RNAseq.pdf", width=9, height=15)
gridExtra::grid.arrange(p1, p2, p3, nrow = 3)
dev.off()

