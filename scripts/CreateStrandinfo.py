import pandas as pd

def CreateStrandinfo(ifile,outdir):
    GFF3 = pd.read_csv(
        filepath_or_buffer=ifile,#'Homo_sapiens.GRCh38.92.ERCC.chr.gtf', 
        sep='\t', 
        header=None,
        names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'],
        skiprows=[i for i in range(25)])

    GFF3 = GFF3[GFF3['source'].notnull()]

    GFF3['gene_names'] = gene_names

    GFF3.index = gene_names

    GFF3['strand'][~GFF3['gene_names'].index.duplicated(keep='first')].to_csv(outdir)#'strandedness.csv')

    #stranded = pd.read_csv('strandedness.csv', header=None, index_col=0)
    print('Wrote Strandedness info to {}'.format(outdir))
    
if __name__ == "__main__":
    import sys
    ifile = sys.argv[1]
#    ifile2 = sys.argv[2]
    outdir = sys.argv[2]
    CreateStrandinfo(ifile,outdir)
