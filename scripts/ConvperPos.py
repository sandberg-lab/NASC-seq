#!/usr/bin/env python
# v0.003

import pysam

def CountConvperPos(bamfile):
    ContigLocs={}#keys= accept_ref_names: values =empty lists }
    AnnoteLocs={}
    for read in bamfile.fetch():
        try:
            if read.get_tag('ST')=='+':
                locs=read.get_tag('TL')
            else:
                locs=read.get_tag('AL')
            if locs[0]!=0:
                if read.reference_name in ContigLocs:
                    ContigLocs[read.reference_name].extend(locs)
                else:
                    ContigLocs[read.reference_name] = list(locs)
                if read.reference_name not in AnnoteLocs:
                    for i,each in enumerate(locs):
                        if i == 0:
                            AnnoteLocs[read.reference_name] = { each :read.get_tag('XT')}
                        else:
                            AnnoteLocs[read.reference_name][each] = read.get_tag('XT')
                else:
                    for i,each in enumerate(locs):
                        if each not in AnnoteLocs[read.reference_name]:
                            AnnoteLocs[read.reference_name][each] = read.get_tag('XT')
        except (ValueError,KeyError):
            continue
    return ContigLocs, AnnoteLocs

def CountReadConverPerConvPos(bam,ContigLocs):
    ConvsPerPos={}
    CoverofPosWithConvs={}
    for key in ContigLocs.keys():
        ContigLocs[key]=sorted(ContigLocs[key])
        ConvsPerPos[key]={}
        k=0
        current=ContigLocs[key][k]
        if len(ContigLocs[key]) == 1:
            ConvsPerPos[key][current] = 1
        else:
            k+=1
            nextone=ContigLocs[key][k]
            while k < len(ContigLocs[key])-1:
                ConvsPerPos[key][current]=1
                while current == nextone and k < len(ContigLocs[key])-1:
                    k+=1
                    nextone=ContigLocs[key][k]
                    ConvsPerPos[key][current]+=1
                current = nextone
                if k < len(ContigLocs[key])-1:
                    k+=1
                    nextone=ContigLocs[key][k]
        CoverofPosWithConvs[key]={}
        for key2 in ConvsPerPos[key].keys():
            try:
                print(key)
                CoverofPosWithConvs[key][key2]=bam.count(key,key2,key2+1)#bam.count(contig=key,start=key2,stop=key2+1)
            except ValueError:
                continue
    return ConvsPerPos,CoverofPosWithConvs

def ExportasVcf(ConvsPerPos,CoverofPosWithConvs, AnnoteLocs):
    #Table Chrom, Pos , ConvsPerPs, CoverofPosWithConvs
    Outputdf =pd.DataFrame(columns=['pos2','convs','covers','chrom','posratio'])
    for key in ConvsPerPos.keys():
        df=pd.DataFrame.from_dict(ConvsPerPos[key], orient='index')#,columns=['pos','convs'])#ConvsPerPos[key])
        df1=pd.DataFrame.from_dict(CoverofPosWithConvs[key], orient='index')#,columns=['pos','covers'])
        df.index.name='pos'
        df1.index.name='pos'
        df.columns = ['convs']
        df1.columns = ['covers']
        df2=df.join(df1)# index='pos')
        df2['pos2'] = df2.index
        df2.index = np.arange(df2.shape[0])
        df2['chrom']=np.repeat(key,df2.shape[0])
        df2['posratio']=df2['convs']/df2['covers']
        df3=pd.DataFrame.from_dict(AnnoteLocs[key], orient='index')
        df3.columns = ['gene_id']
        df2=df2.join(df3, on='pos2')
        Outputdf=Outputdf.append(df2)
    return Outputdf.reset_index(drop=True)

def createTag(d):
    return ''.join([''.join(key) + str(d[key]) + ';' for key in d.keys()])[:-1]

def convInRead(read, qual = 20):
    #vcf_reader = vcf_obj
    specific_conversions = {}
    total_content = {'a' : 0, 'c' : 0, 'g' : 0, 't' : 0}
    specific_conversions[('c', 'A')] = 0
    specific_conversions[('g', 'A')] = 0
    specific_conversions[('t', 'A')] = 0
    specific_conversions[('a', 'C')] = 0
    specific_conversions[('g', 'C')] = 0
    specific_conversions[('t', 'C')] = 0
    specific_conversions[('a', 'G')] = 0
    specific_conversions[('c', 'G')] = 0
    specific_conversions[('t', 'G')] = 0
    specific_conversions[('a', 'T')] = 0
    specific_conversions[('c', 'T')] = 0
    specific_conversions[('g', 'T')] = 0
    specific_conversions[('a', 'N')] = 0
    specific_conversions[('c', 'N')] = 0
    specific_conversions[('g', 'N')] = 0
    specific_conversions[('t', 'N')] = 0

    accept_ref_names= ['chrX','chrY']
    for i in range(1,23):
        accept_ref_names.append('chr'+str(i))
    #if 'chr'+read.reference_name in accept_ref_names:
## has weird addition of 'chr' - format changed in hg38
     #   vcf_iter = vcf_reader.fetch(chrom='chr'+read.reference_name, start=read.reference_start, end=read.reference_end)
     #   snp_list = [record.POS-1 for record in vcf_iter]
   # else:
   #     snp_list=[]
    tC_loc = []
    aG_loc = []

    try:
        refseq = read.get_reference_sequence().lower()
    except (UnicodeDecodeError):
        refseq=''

    for base in total_content.keys():
        total_content[base] += refseq.count(base)
    for pair in read.get_aligned_pairs(with_seq=True):
        try:
            if pair[0] is not None and pair[1] is not None and pair[2] is not None:
                if str(pair[2]).islower() and not read.query_qualities[pair[0]] < qual:
                    specific_conversions[(pair[2],read.seq[pair[0]])] += 1
                    if (pair[2],read.seq[pair[0]]) == ('t', 'C'):
                        tC_loc.append(pair[1])
                    if (pair[2],read.seq[pair[0]]) == ('a', 'G'):
                        aG_loc.append(pair[1])
        except (UnicodeDecodeError, KeyError):
            continue
    SC_tag = createTag(specific_conversions)
    TC_tag = createTag(total_content)
    #print(tC_loc)
    if len(tC_loc) == 0:
        tC_loc.append(0)
    if len(aG_loc) == 0:
        aG_loc.append(0)
    return SC_tag, TC_tag, tC_loc, aG_loc

def addTags(bamfilename, outputname):
    bamfile = pysam.AlignmentFile(bamfilename, 'rb')
    #vcf_reader = vcf.Reader(filename=vcf_filename, compressed=True)
    mod_bamfile = pysam.AlignmentFile(outputname, mode='wb',template=bamfile)
    strandedness = pd.read_csv(strandednessfile, header=None, index_col=0)
    for read in bamfile.fetch():
        try:
            tags = convInRead(read)#, vcf_reader)
            #read.set_tags([('SC',tags[0],'Z'),('TC',tags[1],'Z'), ('TL',tags[2]),('AL',tags[3])])
            read.set_tag('SC',tags[0],'Z')
            read.set_tag('TC',tags[1],'Z')
            read.set_tag('TL',tags[2])#,'B')
            read.set_tag('AL',tags[3])#,'B')
            read.set_tag('ST',strandedness.loc[read.get_tag('XT')][1])
            mod_bamfile.write(read)
        except (ValueError,KeyError):
            continue
    print('Wrote tags to {}'.format(outputname))

if __name__ == "__main__":
    import os
    import pysam
    #import vcf
    import numpy as np
    import pandas as pd
    import sys
    ifile = sys.argv[1]
    outdir = sys.argv[2]
    outdir1 = sys.argv[3]
    cell_id = sys.argv[4]
    strandednessfile = sys.argv[5]
    print('Adding tags to {}'.format(cell_id))
    addTags(ifile,outdir)
    #index outputed bamfile out of pipeline for lazyiness reason

    import os
    import subprocess
    print('Indexing {}'.format(outdir))
    cmd=['samtools index',outdir]
    def run_cmd(cmd):
        subprocess.call(' '.join(cmd),shell=True)
    run_cmd(cmd)

    bam =pysam.AlignmentFile(outdir, 'rb')

    print('Obtaining conversion positions for {}'.format(cell_id))
    ContigLocs, AnnoteLocs=CountConvperPos(bam)

    bam =pysam.AlignmentFile(outdir, 'rb')
    print('Obtaining coverage over conversion position for {}'.format(cell_id))
    ConvsPerPos,CoverofPosWithConvs = CountReadConverPerConvPos(bam,ContigLocs)
    A=ExportasVcf(ConvsPerPos,CoverofPosWithConvs,AnnoteLocs)
    A['cell_id']  = cell_id
    print('Saving result for {} to {}'.format(cell_id, outdir1))
    A.to_csv(outdir1)
