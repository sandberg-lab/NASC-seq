#!/mnt/kauffman/hendgert/Programs/R/R-3.5.1/bin/Rscript

## Alternative to VCF file usage. Instead using conversion rates per position...
## V0.000001
## Usage vcfFilter filedir outdir posratioCutoff cellnumberCutoff

library(data.table)
library(dplyr)
library(tidyr)
library(multidplyr)
library(gplots)

args <- commandArgs(trailingOnly=TRUE)

VCF_merge <- function(filelocations){
    x <- 0
    for (file in Sys.glob(filelocations)){
        print(file)
        if (x==0){
            dt <- fread(file)[,-1]
            x <- x+1
        } else {
            dt <- rbind(dt,fread(file)[,-1])
        }
    }
    return(dt)
}

files <- paste0(args[1],'/*/*_PosTag.csv')
dt <- (VCF_merge(files))
dt$ll <- with(dt, paste(chrom, pos2,sep="_"))
dt <- dt %>% group_by(ll) %>% arrange(desc(posratio)) %>% mutate(pos_index=row_number()) %>% ungroup()
dt <- dt %>% group_by(ll) %>% mutate(cellnumber=length(unique(cell_id))) %>% ungroup()
dt2 <- dt
dt <- dt %>% group_by(ll) %>% mutate(location_sum=sum(posratio)) %>% ungroup()
dt <- select(dt,ll,posratio,pos_index,cellnumber,location_sum)
dt <- spread(dt,key=pos_index,value=posratio)
dt <- dt[order(dt$cellnumber,dt$location_sum,decreasing=TRUE,na.last=TRUE),]
dt <- dt[dt$cellnumber>3,]

outpdf <- paste0(args[2],'/ConversionPerPosition_Top2000.pdf')
pdf(outpdf,height=10,width=6)
heatmap.2(as.matrix(dt[1:2000,-c(1:3)]),Colv=NA,Rowv=NA,trace='none',ylab='Position (sorted by detection and sum of conversion rates)',xlab='Cell (sorted by conversion ratios)',dendrogram='none',labRow='',labCol='',key.title = 'Conversion rate',main = 'Conversions on top 2000 positions')
dev.off()
pos.sel <- unique(dt2[dt2$posratio >= as.numeric(args[3]) & dt2$cellnumber > as.numeric(args[4]),c('chrom','pos2')])
outcsv <- paste0(args[2],'/posfile.csv')
write.csv(pos.sel,outcsv)