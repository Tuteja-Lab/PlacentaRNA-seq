# Install page:

# CRAN PACKAGES


# BIOCCONDUCTOR PACKAGES
install.packages(c('dplyr', 'tidyverse', 'ggplot2', 'BiocManager', 'dendextend', 'matrixStats',
                   'devtools', 'pca3d', 'scatterplot3d', 'pheatmap', 'dendextend', 'remotes'))


require(devtools)
install_version("stats", version = "3.6.3", repos = "http://cran.us.r-project.org")

BiocManager::install(c('biomaRt' ))
BiocManager::install('stats', version = '3.6.3')
BiocManager::install('sleuth', version = '0.30.0')
BiocManager::install('org.Mm.eg.db', version = '3.13.0')


library("kohonen")
library("RclusTool")
library("tximport")
library("clusterProfiler")
library("org.Mm.eg.db")
library("GENIE3", 1.16)

# Versions:
#   hclust, prcomp:  package stats [114], version 3.6.3
# phyper() (package stats, version 3.6.3
#           kohonen [115], version 3.0.10
#           RclusTool [116], version 0.91.3
#           STRING database (version 11.0
#                            tximport (version 1.14.2
#                                      GENIE3 (version 1.8.0) or 1.16
#                                      ClusterProfiler (version 4.0.5
#                                                       
#                                                       LinSeed (version 0.99.2
# 
?install_github
# Linseed
# install.packages("remotes")
# remotes::install_github("ctlab/LinSeed")