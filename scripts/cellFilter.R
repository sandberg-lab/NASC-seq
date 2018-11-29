# V0.002
## Usage cellFilter.R infile outfile totalReads assignedReads percentageAssigned numberofHighExprGenes
## Produces list of cells that pass QC filter

args <- commandArgs(trailingOnly=TRUE)
qcdata <- data.frame(readRDS(args[1]))
numBefore <- nrow(qcdata)
qcdata <- qcdata[as.numeric(qcdata[,'totalReads']) >= as.numeric(args[3]),]
qcdata <- qcdata[as.numeric(qcdata$Assigned) >= as.numeric(args[4]),]
qcdata <- qcdata[as.numeric(qcdata$percentageAssigned) >= as.numeric(args[5]),]
numAfter <- nrow(qcdata)
cellsPassing <- rownames(qcdata)
write(cellsPassing,file=args[2],ncolumns=1,append=FALSE)

print(paste(numAfter,'out of',numBefore, 'cells passed filtering...'))
