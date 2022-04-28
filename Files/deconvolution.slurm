#!/bin/bash

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=128G # Maximum amount of memory this job will be given, try to estimate this to the best of your ability. This asks for 128 GB of ram.
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mail-user=krkies@iastate.edu  # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="/work/LAS/geetu-lab/krkies/PlacentaRNA-seq/deconvolution_422.out" 
#SBATCH --error="/work/LAS/geetu-lab/krkies/PlacentaRNA-seq/deconvolution_422.errors" 



module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq
module load gsl/2.5-fpqcpxf
module load udunits/2.2.24-yldmp4h
module load gdal/2.4.4-nw2drgf
module load geos/3.8.1-2m7gav4

export R_LIBS_USER=/work/LAS/geetu-lab/hhvu/R/x86_64-pc-linux-gnu-library/4.0
mkdir -p $R_LIBS_USER


Rscript deconvolutionAnalysis.R
