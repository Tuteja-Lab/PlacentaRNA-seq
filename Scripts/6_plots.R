library(ggplot2)

dir <- "Z:/hhvu/Project1_2/RNA-seq/6_clusterProfiler/GENIE3/"
permResults <- "perm-max1000/"
rObjects <- "max1000/"
timepoint <- "e9.5"
subnet <- "e9.5_3"
load(paste0(dir, rObjects, subnet, "_BP_filtered.rda"))

original <- goTerms
#interestingTerms <- original[grep("proliferation|immune|immunity|cell migration|cell motility|branching|labyrinth|embryo|vascular|vasculature|angiogenesis|blood vessel|trophoblast|tube|hypoxia|placenta|vitamin transport|gas transport|nutrient", go2$Description),]

files <- list.files(paste0(dir, permResults))
files <- files[grep(subnet, files)]
files <- files[grep("txt", files)]
terms <- sapply(strsplit(files, "_"), function(x) x[3])
terms <- terms[-1]

plots <- c()
for (i in 1:length(terms)) {
  term <- goTerms[goTerms$Description == terms[i],]
  if (length(files[grep(terms[i], files)]) > 1) {
    temp <- files[grep(terms[i], files)][1]
    permRes <- read.table(paste0(dir, permResults, temp), header = T)
  }  else {
    permRes <- read.table(paste0(dir, permResults, files[grep(terms[i], files)]), header = T)
  }
  colnames(permRes) <- "permRes"
  p <- ggplot(permRes, aes(x = log10(permRes))) + geom_histogram() +
    geom_vline(data = term, aes(xintercept = -log10(qvalue)), color="blue", linetype="dashed", size=1) +
    labs(title = paste0("\'", term$Description, "\'", "\n"), x = "-log10(qvalue)", y = "Count") +
    theme(plot.title = element_text(size=8))
  plots[[i]] <- p
}
#, "in ", subnet, " network"


pdf(paste0(dir, permResults, subnet, "_", "permutationResults2.pdf"), width = 15, height = 8)
gridExtra::grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
                        plots[[5]], plots[[6]], plots[[7]], plots[[8]],
                        plots[[9]], plots[[10]], plots[[11]], #plots[[12]],
                        #plots[[13]], #plots[[14]], plots[[15]], plots[[16]],
                        #plots[[17]], plots[[18]], plots[[19]], plots[[20]],
                        #plots[[21]], plots[[22]], plots[[23]], #plots[[24]],
                        #plots[[25]], plots[[26]], plots[[27]], plots[[28]],
                        #plots[[29]], plots[[30]], plots[[31]],
                        nrow = 4, top = grid::textGrob(paste0("Subnetwork ", subnet, " STRING"), gp=grid::gpar(fontsize=20,font=3)))
dev.off()
