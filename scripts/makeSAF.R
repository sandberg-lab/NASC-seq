#!/mnt/kauffman/hendgert/Programs/R/R-3.5.1/bin/Rscript

## Intron annotation feature from ZUMIs

library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(GenomicFeatures)
library(GenomicRanges)

args <- commandArgs(trailingOnly=TRUE)


makeSAF<-function(gtf){
  print("Loading reference annotation from:")
  print(gtf)
  txdb <- suppressWarnings(suppressMessages(GenomicFeatures::makeTxDbFromGFF(file=gtf, format="gtf")))
  
  ## Make Gene-range GR-object
  se <- suppressMessages(
    AnnotationDbi::select(txdb, keys(txdb, "GENEID"),
                          columns=c("GENEID","TXCHROM","TXSTART","TXEND","TXSTRAND"),
                          keytype="GENEID") %>%
      dplyr::group_by(GENEID,TXCHROM,TXSTRAND)  %>% 
      dplyr::mutate( txstart =ifelse(TXSTART<TXEND,min(TXSTART),min(TXEND)),
                     txend  =ifelse(TXSTART<TXEND,max(TXEND),min(TXSTART) ) ) %>%
      dplyr::select(GENEID,TXCHROM,TXSTRAND,txstart,txend)  %>% unique()
  )
  
  gr.gene<-GenomicRanges::GRanges(seqnames = se$TXCHROM,
                                  ranges =  IRanges(start= se$txstart,
                                                    end=  se$txend,
                                                    names=se$GENEID),
                                  strand =  se$TXSTRAND,
                                  gid    =  se$GENEID)
  
  ### Get non-overlapping Introns/Exons
  intron<-GenomicFeatures::intronsByTranscript(txdb, use.names=T)
  exon<-GenomicFeatures::exonsBy(txdb, by="tx",use.names=T)
  
  intron.exon.red <- c( GenomicRanges::reduce(unlist(intron),ignore.strand=T), GenomicRanges::reduce(unlist(exon),ignore.strand=T) )
  intron.exon.dis <- GenomicRanges::disjoin(intron.exon.red, ignore.strand=T)
  intron.only<-GenomicRanges::setdiff(intron.exon.dis, unlist(exon) ,ignore.strand=T)
  
  ol.in<-GenomicRanges::findOverlaps(intron.only, gr.gene, select="arbitrary")
  ol.ex<-GenomicRanges::findOverlaps(unlist(exon), gr.gene, select="arbitrary")
  
  intron.saf<-data.frame(GeneID= names(gr.gene)[ol.in],
                         Chr   = seqnames(intron.only),
                         Start = start(intron.only),
                         End	 =   end(intron.only),stringsAsFactors = F)
  exon.saf<-data.frame(GeneID= names(gr.gene)[ol.ex],
                       Chr   = seqnames(unlist(exon)),
                       Start = start(unlist(exon)),
                       End	 =   end(unlist(exon)),
                       Strand =  strand(unlist(exon)),stringsAsFactors = F)
  
  intron.saf<-dplyr::left_join(intron.saf,unique(exon.saf[,c("GeneID","Strand")]),by=c("GeneID"))
  
  saf <- list(introns=intron.saf,exons=exon.saf)
  print("Annotation loaded!")
  #  safout <- paste(out,"/zUMIs_output/expression/",sn,".annotationsSAF.rds",sep="")
  #  saveRDS(saf, file=safout)
  rm(se,gr.gene,intron,exon,intron.exon.red,intron.exon.dis,intron.only,ol.ex,ol.in,intron.saf,exon.saf)
  return(saf)
}


fileIn <- args[1]
fileOut <- args[2]
intron.saf <- makeSAF(fileIn)
saveRDS(intron.saf,fileOut)