#!/usr/bin/env python
##

# G.J. Hendriks

## V0.1000
pipelineVersion=0.1000

import os, argparse, sys, subprocess, time, glob, ntpath
from joblib import Parallel, delayed
import pandas as pd
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experimentdir', required=True)
    parser.add_argument('-p', '--numCPU', required=True, type=int)
    parser.add_argument('-f', '--flag',type=str)
    o = parser.parse_args()

def safe_mkdir(f):
  if not os.path.exists(f):
    os.mkdir(f)

def run_cmd(cmd):
    subprocess.call(" ".join(cmd),shell=True)

def cellListCheck(cellname,experimentdir):
	text_file=open(os.path.join(experimentdir,'QC/cellFilter/cellFilter.txt'))
	lines = text_file.read().split('\n')
	return cellname in lines

def cellNameExtration(filename):
	cellname='_'.join(os.path.basename(cell).split('_')[:3])
	return cellname

def makePickleList(directory,file):
        pklFileList = glob.glob(os.path.join(directory,'*.pkl'))
        with open(file, 'w') as f:
            for pkl in pklFileList:
                f.write("%s\n" % pkl)

execfile(os.path.join(o.experimentdir,'config.py'))
print ('config.py found in: %r')%os.path.join(o.experimentdir,'config.py')

## Genome versions and settings
bases = experimentInfo['readlength']
verbose = experimentInfo['verbose']
gtf = experimentInfo['gtf']
readcountdir=experimentInfo['readcountdir']
output=experimentInfo['output']
aligndir=experimentInfo['aligndir']
gnv=experimentInfo['gnv']
gtf=experimentInfo['gtf']
sjdb=experimentInfo['sjdb']
strandFile=experimentInfo['strandFile']
stanFile=experimentInfo['stanFile']

## Distributions of dependencies
rootDir=distributions['rootDir']
trimgaloreDist=distributions['trimgaloreDist']
starDist=distributions['starDist']
picardDist=distributions['picardDist']
javaDist=distributions['javaDist']
RscriptDist=distributions['RscriptDist']
## Standard directory/filenames
fastqfiles='fastqFiles'
trimfiles='trimmedFastqFiles'
bamfiles='bamFiles'
aligndir='aligned_bam'
removedup='duplRemoved_bam'
annotatedfiles='annotated_bam'
intronannotatefiles='intronAnnotated_bam'
annotatedsortedfiles='annotated_sorted_bam'
qcfiles='QC'
positionfiles='position_counts'
taggedfiles='tagged_bam'
featurecountQC='featureCount_QC'
cellFilter='cellFilter'
vcfFilter='vcfFilter'
filterTaggedFiles='filteredTagged_bam'
p_e='p_e'
errorRates='errorRates'
outfiles='outfiles'
pi_g='pi_g'
pi_g_specific='pi_g_specific'
pkl_files = 'pkl_files'

## Commandlogfile
commandlogfile = open(os.path.join(o.experimentdir,'commandlog.txt'), 'a')
commandlogfile.write('PipelineVersion: %f\n' % pipelineVersion)

print ('read length is %r') %bases
print ('verbose is %r')%verbose

if o.flag=='trim' or o.flag=='all':
	print ('Performing trimming using TrimGalore-0.4.5...')
	safe_mkdir(os.path.join(o.experimentdir,fastqfiles,trimfiles))
	cmds = []
	for cell in glob.glob(os.path.join(o.experimentdir,fastqfiles,'rawdata','*')):
		cellname=os.path.basename(cell)
		files = glob.glob(os.path.join(cell,'*'))
		outdir = os.path.join(cell,'..','..',trimfiles,cellname)
		safe_mkdir(outdir)
		cmd = [trimgaloreDist, '--nextera --length 20', files[0],files[1],'-o',outdir, '--paired --dont_gzip']
		cmds.append(cmd)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)

if o.flag=='align' or o.flag=='all':
  safe_mkdir(os.path.join(o.experimentdir,bamfiles))
  safe_mkdir(os.path.join(o.experimentdir,bamfiles,aligndir))
  cmd1 = [starDist,'--genomeLoad LoadAndExit','--genomeDir',gnv,'--limitBAMsortRAM 150000000000']
  cmds = []
  for cell in glob.glob(os.path.join(o.experimentdir,fastqfiles,trimfiles,'*')):
    cellname='_'.join(os.path.basename(cell).split("_")[:3])
    safe_mkdir(os.path.join(o.experimentdir,bamfiles,aligndir,cellname))
    files = glob.glob(os.path.join(cell,'*.fq'))
    outprefix = os.path.join(o.experimentdir,bamfiles,aligndir,cellname,cellname+"_")
    if os.path.isfile(outprefix+'Aligned.sortedByCoord.out.bam'):
      print(outprefix+'Aligned.sortedByCoord.out.bam already exists... Not performing alignment for this cell')
    else:
      cmd = [starDist,'--runThreadN 4','--outFileNamePrefix',outprefix,'--genomeDir',gnv,
      '--readFilesIn',files[0],files[1],'--alignSJoverhangMin 1000 --alignSJDBoverhangMin 1 --bamRemoveDuplicatesType UniqueIdentical --outFilterMismatchNoverReadLmax 1',
      '  --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 --outSAMattributes MD --outSAMtype BAM SortedByCoordinate',
      ' --scoreDelOpen -10000 --scoreInsOpen -10000 --outSAMmultNmax 1 --genomeLoad LoadAndKeep --limitBAMsortRAM 150000000000']
      cmds.append(cmd)
  run_cmd(cmd1)
  Parallel(n_jobs=int(o.numCPU/4))(delayed(run_cmd)(cmd) for cmd in cmds)
  for cmd in cmds:
    commandlogfile.write('%s\n' % cmd)
  print ("All fastq files have been aligned to genome found in %r") %gnv

if o.flag=='removegenome' or o.flag=='all':
  cmd = [starDist,'--genomeLoad Remove','--genomeDir', gnv]
  run_cmd(cmd)
  print ("Succesfully removed genome from shared memory")

if o.flag=='removeduplicates' or o.flag=='all':
	print ('Removing duplicates using Picard MarkDuplicates function...')
	cmds = []
	outdir=os.path.join(o.experimentdir,bamfiles,removedup)
	safe_mkdir(outdir)
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,aligndir,"*/*sortedByCoord.out.bam")):
		cellname='_'.join(os.path.basename(cell).split("_")[:3])
		safe_mkdir(os.path.join(outdir,cellname))
		outfile=os.path.join(outdir,cellname,cellname+"_removeDupl.bam")
		metricsfile=os.path.join(outdir,cellname,cellname+"_removeDuplMetrics.txt")
		cmd = [javaDist,' -XX:ParallelGCThreads=2 -jar', picardDist ,'  MarkDuplicates', 'I='+cell, 'O='+outfile,'M='+metricsfile, 'REMOVE_DUPLICATES=true']
		cmds.append(cmd)
	if verbose: print(" ".join(cmds))
	Parallel(n_jobs=int(o.numCPU/2))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Duplicate reads have been removed")

if o.flag=='removeduplicates' or o.flag=='all':
  cmds = []
  for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,removedup,'*/*_removeDupl.bam')):
    cmd = ['samtools index',cell]
    cmds.append(cmd)
  if verbose: print(" ".join(cmds))
  Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
  for cmd in cmds:
    commandlogfile.write('%s\n' % cmd)
  print ("Indexing of bam files has been completed")

if o.flag=='annotate' or o.flag=='all':
	cmds = []
	print ('Annotating features using FeatureCounts from RSubRead...')
	safe_mkdir(os.path.join(o.experimentdir,bamfiles,annotatedfiles))
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,removedup,'*/*_removeDupl.bam')):
		cellname='_'.join(os.path.basename(cell).split("_")[:3])
		outpath = os.path.join(o.experimentdir,bamfiles,annotatedfiles,cellname)
		safe_mkdir(outpath)
		cmd = [RscriptDist,os.path.join(rootDir,'scripts/RsubReadFeatureCounts.R'),cell, gtf, outpath]
		cmds.append(cmd)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("First pass Feature counting and BAM annotation finished")

if o.flag=='annotate' or o.flag=='all':
	print('Sorting and indexing using Samtools...')
	safe_mkdir(os.path.join(o.experimentdir,bamfiles,annotatedsortedfiles))
	cmds = []
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,annotatedfiles,'*/*featureCounts.bam')):
		cellname='_'.join(os.path.basename(cell).split("_")[:3])
		safe_mkdir(os.path.join(o.experimentdir,bamfiles,annotatedsortedfiles,cellname))
		outfile=os.path.join(o.experimentdir,bamfiles,annotatedsortedfiles,cellname,cellname+'_sorted.bam')
		cmd = ['samtools sort',cell,'-o',outfile]
		cmds.append(cmd)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	cmds = []
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,annotatedsortedfiles,'*/*sorted.bam')):
		cmd = ['samtools index',cell]
		cmds.append(cmd)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Indexing of bam files has been completed")

if o.flag=='conversiontag' or o.flag=='all':
	print ('Adding conversion information to tags in bamfiles')
	cmds = []
	safe_mkdir(os.path.join(o.experimentdir,bamfiles,taggedfiles))
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,annotatedsortedfiles,'*/*sorted.bam')):
		cellname = '_'.join(os.path.basename(cell).split('_')[:3])
		safe_mkdir(os.path.join(o.experimentdir,bamfiles,taggedfiles,cellname))
		outfile = os.path.join(o.experimentdir,bamfiles,taggedfiles,cellname,cellname+'_PositionTagged.bam')
		outfile2 = os.path.join(o.experimentdir,bamfiles,taggedfiles,cellname,cellname+'_PosTag.csv')
		cmd = ['python2.7', os.path.join(rootDir,'scripts/ConvperPos.py'),cell,outfile,outfile2,cellname,strandFile]
		cmds.append(cmd)
	if verbose: print(" ".join(cmds))
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Conversions have been added to bamfile tags, and position information has been stored")

if o.flag=='vcfFilter':
	print('Selecting positions that are possible SNPs based on high positional conversion rates...')
	safe_mkdir(os.path.join(o.experimentdir,qcfiles))
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,vcfFilter))
	indir = os.path.join(o.experimentdir,bamfiles,taggedfiles)
	outdir = os.path.join(o.experimentdir,qcfiles,vcfFilter)
	posratioCutoff = cutoffs['posratioCutoff']
	cellnumberCutoff = cutoffs['cellnumberCutoff']
	cmd=[RscriptDist,os.path.join(rootDir,'scripts/vcfFilter.R'),indir,outdir,str(posratioCutoff),str(cellnumberCutoff)]
	commandlogfile.write('%s\n' % cmd)
	run_cmd(cmd)

if o.flag=='tagFilter':
	cmds=[]
	print ('Refining positions considered as conversions based on snp-filter and paired-end overlaps...')
	safe_mkdir(os.path.join(o.experimentdir,bamfiles,filterTaggedFiles))
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,taggedfiles,'*/*_PositionTagged.bam')):
		cellname='_'.join(os.path.basename(cell).split('_')[:3])
		safe_mkdir(os.path.join(o.experimentdir,bamfiles,filterTaggedFiles,cellname))
		outfile=os.path.join(o.experimentdir,bamfiles,filterTaggedFiles,cellname,cellname+'_taggedFiltered.bam')
		inposfile=os.path.join(o.experimentdir,qcfiles,vcfFilter,'posfile.csv')
		cmd=['python2.7', os.path.join(rootDir,'scripts/filter_reads_paired.py'),cell,outfile,inposfile]
		cmds.append(cmd)
	if verbose: print(" ".join(cmds))
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Filtered conversions to remove potential SNPs and overcounting")

if o.flag=='tagFilter':
	cmds=[]
	print ("Indexing filtered and tagged bamfiles...")
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,filterTaggedFiles,'*/*_taggedFiltered.bam')):
		cmd = ['samtools index',cell]
		cmds.append(cmd)	
	if verbose: print(" ".join(cmds))
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Indexing finished")

if o.flag=='cellQC':
	print ('Performing quality control for all cells based on featureCounts QC data...')
	cmds=[]
	safe_mkdir(os.path.join(o.experimentdir,qcfiles))
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,featurecountQC))
	rdslocationpath=os.path.join(o.experimentdir,bamfiles,annotatedfiles)
	outdir=os.path.join(o.experimentdir,qcfiles,featurecountQC)
	print ('rds location path = %r')%rdslocationpath
	print ('outdir = %r')%outdir
	cmd=[RscriptDist,os.path.join(rootDir,'scripts/cellQC.R'),rdslocationpath,outdir]
	commandlogfile.write('%s\n' % cmd)
	run_cmd(cmd)

if o.flag=='cellFilter':
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,))
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,cellFilter))
	print ('Filtering cells based on filter variables in config.py')	
	cmds=[]
	infile=os.path.join(o.experimentdir,qcfiles,featurecountQC,'QC_data_featureCounts.rds')
	outfile=os.path.join(o.experimentdir,qcfiles,cellFilter,'cellFilter.txt')
	totalReads=cutoffs['totalReads']
	assignedReads=cutoffs['assignedReads']
	percentageAssigned=cutoffs['percentageAssigned']
	cmd=[RscriptDist,os.path.join(rootDir,'scripts/cellFilter.R'),infile,outfile,str(totalReads),str(assignedReads),str(percentageAssigned)]
	commandlogfile.write('%s\n' % cmd)
	run_cmd(cmd)

if o.flag=='calculatePE':
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,p_e))
	safe_mkdir(os.path.join(o.experimentdir,qcfiles,errorRates))
	print ('Calculating PE for all cells that passed filter...')
	cmds=[]
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,filterTaggedFiles,'*/*_taggedFiltered.bam')):
		cellname=cellNameExtration(cell)
		if cellListCheck(cellname,o.experimentdir):
			outfile=os.path.join(o.experimentdir,qcfiles,errorRates,cellname+'_ErrorRates.csv')
			outfile2=os.path.join(o.experimentdir,qcfiles,p_e,cellname+'_p_e.txt')
			cmd=['python2.7',os.path.join(rootDir,'scripts/PECalculation_GJH.py'),cell,outfile,outfile2]
			cmds.append(cmd)
	if verbose: print(" ".join(cmds))
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Pe has been calculated for all cells that passed filter.")

if o.flag=='prepareData':
	pickle_outdir = os.path.join(o.experimentdir,outfiles,pkl_files)
	safe_mkdir(os.path.join(o.experimentdir,outfiles))
	safe_mkdir(pickle_outdir)
	cmds = []
	for cell in glob.glob(os.path.join(o.experimentdir,bamfiles,filterTaggedFiles,'*/*_taggedFiltered.bam')):
		cellname = cellNameExtration(cell)
		if cellListCheck(cellname,o.experimentdir):
			pe_file = open(os.path.join(o.experimentdir,qcfiles,p_e,cellname+'_p_e.txt'))
			pe_estimate = pe_file.read().split('\n')[0]
			outfile = os.path.join(pickle_outdir,cellname+'_prepared.pkl')
			cmd = ['python3',os.path.join(rootDir,'scripts/prepare_pickles.py'),cell,outfile,pe_estimate,cellname]
			cmds.append(cmd)
	print (cmds)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
	for cmd in cmds:
		commandlogfile.write('%s\n' % cmd)
	print ("Data has been prepared for scalable processing...")

if o.flag=='processData':
	safe_mkdir(os.path.join(o.experimentdir,outfiles,'outPickles'))
	pklfile = os.path.join(o.experimentdir,outfiles,'PklList.txt')
	indir = os.path.join(o.experimentdir,outfiles,'pkl_files')
	makePickleList(indir,pklfile)
	outfile = os.path.join(o.experimentdir,outfiles,'outPickles/pi_g_results.pkl')
	logfile_pi_g = os.path.join(o.experimentdir,outfiles,'outPickles/logfile.txt')
	cmd=['python3', os.path.join(rootDir,'scripts/estimate_pi_g_STAN_forAWS_v2.py'),outfile,pklfile,stanFile,int(o.numCPU),'>',logfile_pi_g]
	run_cmd(cmd)

if o.flag=='summarize':
	safe_mkdir(os.path.join(o.experimentdir,outfiles,'summary'))
	readcountRDS = os.path.join(o.experimentdir,qcfiles,featurecountQC,'Counttable*.rds')
	readcountOut = os.path.join(o.experimentdir,outfiles,'readCounts.csv')
	cmd0 = [RscriptDist, os.path.join(rootDir,'scripts/convertRdsCountable.R'),readcountRDS,readcountOut]
	run_cmd(cmd0)
	indir = os.path.join(o.experimentdir,outfiles,'outPickles','*')
	outfile = os.path.join(o.experimentdir,outfiles,'resultfile_pickles.txt')
	makePickleList(indir,outfile)
	outfile2 = os.path.join(o.experimentdir,outfiles,'summary')
	readcounts = os.path.join(o.experimentdir)
	cmds=[]
	cmd1 = ['python3', os.path.join(rootDir,'scripts/separateTranscriptomes_dicts.py'),outfile,readcountOut,outfile2]
	cmds.append(cmd1)
	cmd2 = ['python3', os.path.join(rootDir,'scripts/getParamTable.py'),outfile,'pi_g','mean',outfile2]
	cmds.append(cmd2)
	cmd3 = ['python3', os.path.join(rootDir,'scripts/getParamTable.py'),outfile,'pi_g','sd',outfile2]
	cmds.append(cmd3)
	cmd4 = ['python3', os.path.join(rootDir,'scripts/getParamTable.py'),outfile,'pi_g','2.5%',outfile2]
	cmds.append(cmd4)
	cmd5 = ['python3', os.path.join(rootDir,'scripts/getParamTable.py'),outfile,'pi_g','97.5%',outfile2]
	cmds.append(cmd5)
	cmd6 = ['python3', os.path.join(rootDir,'scripts/getModeTable.py'),outfile,outfile2]
	cmds.append(cmd6)
	Parallel(n_jobs=int(o.numCPU))(delayed(run_cmd)(cmd) for cmd in cmds)
