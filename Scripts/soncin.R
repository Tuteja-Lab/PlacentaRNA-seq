library("dplyr")
library("lumi")
library("ggplot2")

id <- read.table("/work/LAS/geetu-lab/hhvu/soncin.et.al.2018/GPL10558-50081.txt", header = T, fill = T, sep = "\t", quote = "")
idmm10 <- read.table("/work/LAS/geetu-lab/hhvu/soncin.et.al.2018/GPL6885-11608.txt", header = T, fill = T, sep = "\t", quote = "")

soncin <- read.table("/work/LAS/geetu-lab/hhvu/soncin.et.al.2018/GSE100051_non-normalized.txt", header = T, sep = "\t")
soncin2 <- soncin[,1:(54*2+1)]

avg <- soncin2[,c(1, grep("SAMPLE", names(soncin2)))]
avg <- inner_join(avg, id[, c("ID", "Species", "ILMN_Gene")], by = c("X" = "ID"))
avg1 <- avg[which(avg$ILMN_Gene %in% avg$ILMN_Gene[duplicated(avg$ILMN_Gene)]),]
avg1[,2:(ncol(avg1)-2)] <- sapply(avg1[,2:(ncol(avg1)-2)], as.numeric)
avg1$avg <- rowMeans(avg1[,2:(ncol(avg1)-2)])

p <- soncin2[,c(grep("X", names(soncin2)))]
p <- inner_join(p, id[, c("ID", "Species", "ILMN_Gene")], by = c("X" = "ID"))
p1 <- p[p$X %in% avg1$X,]
p1[,2:(ncol(p1)-2)] <- sapply(p1[,2:(ncol(p1)-2)], as.numeric)
p1$avg <- rowMeans(p1[,2:(ncol(p1)-2)])

avg1$p.avg <- p1$avg

keep <- data.frame(matrix(data=NA, ncol = ncol(avg1)))
names(keep) <- names(avg1)
for (i in unique(avg1$ILMN_Gene)) {
  sub <- avg1[avg1$ILMN_Gene == i,]
  sub <- sub[intersect(which(sub$avg == max(sub$avg)), which(sub$p.avg == min(sub$p.avg))),]
  keep <- rbind(keep, sub)
}
keep <- keep[!is.na(keep$X),]

p2 <- p[p$ILMN_Gene %in% setdiff(p$ILMN_Gene, p$ILMN_Gene[duplicated(p$ILMN_Gene)]),]
p2 <- rbind(p2, p[p$X %in% keep$X,])
p2[,2:(ncol(p2)-2)] <- sapply(p2[,2:(ncol(p2)-2)], as.numeric)
p2$p.avg <- rowMeans(p2[,2:(ncol(p2)-2)])
p2 <- p2[p2$p.avg < 0.01,]

avg2 <- avg[avg$X %in% p2$X,]

rownames(avg2) <- avg2$ILMN_Gene
avg2 <- avg2[,2:55]
avg2 <- as.matrix(avg2)
avg3 <- rsn(avg2)

plot(density(log2(avg3[,1])),col="blue",lwd=3)

a <- reshape2::melt(avg3)
ggplot(aes(x=log2(value), colour=Var2), data=a) + geom_density()

ggplot(aes(x=log2(value), colour=Var2), data=a[a$Var2 %in% c("X.1", "X.2", "X.3"),]) + geom_density() + ggtitle("PLACENTA_Wk4_Tr1")

tpm2 <- read.table("/work/LAS/geetu-lab/hhvu/project1_2/rna-seq/PlacentaRNA-seq/Files/tpmForClustering.txt", header = T)
t2g <- read.table("Files/t2g.txt", header = T)

load(file = "Files/combine-test-expression1.Rdata")
mouseHumanOrthologs<-dataset$GRCH38$mouseHumanOrthologs
hvGenes <- t2g[t2g$target_id %in% rownames(tpm2),]
hvGenes2 <- mouseHumanOrthologs[mouseHumanOrthologs$Gene.stable.ID %in% hvGenes$ens_gene,]

length(intersect(hvGenes2$Human.gene.name, avg2$ILMN_Gene))
intersect(hvGenes2$Human.gene.name, avg2$ILMN_Gene)

length(unique(hvGenes$ext_gene))
length(unique(hvGenes2$Gene.name))
