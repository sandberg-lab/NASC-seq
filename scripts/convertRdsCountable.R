#!/mnt/kauffman/hendgert/Programs/R/R-3.5.1/bin/Rscript
# V0.001
## Usage cellFilter.R infile outfile totalReads assignedReads percentageAssigned numberofHighExprGenes
## Produces list of cells that pass QC filter

args <- commandArgs(trailingOnly=TRUE)
readcounts <- readRDS(args[1])
write.csv(t(readcounts),args[2])
print(paste('Wrote RDS file',args[1],'to file',args[2]))
