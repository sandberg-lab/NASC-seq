
# This resolves issue/enhancement #4 from Sandberg lab NASC-seq repository
# 
# usage: rscript ./strandednessFilePreparation.R gtf-file strandedness-file
# 
# goal was to have a function that prepares the strandedness file.
# 
# mouse:
# example gtf: /home/hendgert/data/20191024/NASC-seqV2_optimization.final_annot.gtf
# example strandedness:/home/hendgert/resources/genomes/mm10_ercc/strandedness_diySpike.csv
# 
# human:
# example gtf: /home/hendgert/resources/genomes/hg38/Homo_sapiens.GRCh38.97.gtf
# example strandedness: /home/hendgert/programs/NASC-seq/data/strandedness.csv
# example usage: /home/hendgert/programs/R-3.5.2/bin/Rscript /home/hendgert/scripts/NASC-seq/NASC-seq/pipeline_11oct18/4sU/NASCseqV2/scripts/strandednessFilePreparation.R /home/hendgert/resources/genomes/hg38/Homo_sapiens.GRCh38.97.gtf /home/hendgert/tempStrandedness.csv


library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

gtf <- fread(args[1])

strandedness <- gtf %>%
  dplyr::select(V7,V9) %>%
  mutate(strand=V7) %>%
  separate(col=V9,into='geneIDpart1',sep=';') %>%
  distinct() %>%
  separate(col=geneIDpart1,into=c(NA,'geneIDpart2'),sep=" ") %>%
  separate(col=geneIDpart2,into=c(NA,'geneID')) %>%
  dplyr::select(geneID,strand)

fwrite(strandedness,args[2],sep=',',col.names = FALSE)
