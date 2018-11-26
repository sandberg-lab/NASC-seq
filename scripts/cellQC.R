#!/mnt/kauffman/hendgert/Programs/R/R-3.5.1/bin/Rscript

#v0.0005

## Usage cellQC rds_file_locations outdir
## Note that rds_file_locations is a single name with wildcards to refer to all rds files!

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)

featureCountsQC <- function(filelocations){
    x <- 0
    for (file in Sys.glob(filelocations)){
        if (x==0){
            print(file)
            obj <- readRDS(file)
            df <- obj[[1]]
            qc <- obj[[4]]
            colnames(df)[1] <- strsplit(file,'/')[[1]][length(strsplit(file,'/')[[1]])-1]## -2
            colnames(qc)[2] <- strsplit(file,'/')[[1]][length(strsplit(file,'/')[[1]])-1]## -2
            x <- x+1
        } else {
            print(file)
            obj <- readRDS(file)
            df2 <- obj[[1]]
            qc2 <- obj[[4]]
            cellname <- strsplit(file,'/')[[1]][length(strsplit(file,'/')[[1]])-1]## -2
            colnames(df2)[1] <- cellname
            colnames(qc2)[2] <- cellname
            df <- cbind(df,df2)
            qc <- cbind(qc,qc2[,cellname])
            colnames(qc)[length(colnames(qc))] <- cellname
            x <- x+1
        }
    }
    return(list(df,qc))
}

qcdir <- paste0(args[1],'/*/*.rds')
print(qcdir)
qclist <- featureCountsQC(qcdir)

counts <- qclist[[1]]
qcdata <- as.data.frame(t(data.frame(qclist[[2]],row.names=1)))
qcdata$totalReads <- apply(qcdata,1,function(x) sum(as.numeric(x)))
qcdata$percentageAssigned <- apply(qcdata,1,function(x) 100*x['Assigned']/x['totalReads'])
qcdata$cutoffPass100 <- apply(counts,2,function(x) sum(x>=100))

p1 <- ggplot(qcdata,aes(x=totalReads)) + geom_histogram() + xlab('Total reads / cell')
p2 <- ggplot(qcdata,aes(x=Assigned)) + geom_histogram() + xlab('Assigned reads / cell')
p3 <- ggplot(qcdata,aes(x=percentageAssigned)) + geom_histogram() + xlab('Percentage of reads aligning / cell')
p4 <- ggplot(qcdata,aes(x=cutoffPass100)) + geom_histogram() +xlab('Number of genes with > 100 reads / cell')
print(paste('argument 1 =',args[1]))
print(paste('argument 2 =',args[2]))
p <- plot_grid(p1,p2,p3,p4,ncol=2,labels=c('A','B','C','D'))
save_plot(paste0(args[2],'/cellQC.pdf'),p)
saveRDS(qcdata,paste0(args[2],'/QC_data_featureCounts.rds'))
saveRDS(counts,paste0(args[2],'/Counttable_featureCounts.rds'))

