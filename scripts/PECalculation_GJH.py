## GJH


## V0.2
## Quick PE calculation script that follows standard data format structure
## Can also be used for STAR optimization visualization

import pandas as pd
import numpy as np
import pysam

## From Anton
def parseSCTag(read):
	tag = read.get_tag('SC')
	splittag = tag.split(';')
	specific_conversions = {}
	for c in splittag:
		specific_conversions[(c[0:2])] = np.int_(c[2:])
	return specific_conversions

def parseTCTag(read):
	tag = read.get_tag('TC')
	splittag = tag.split(';')
	total_content = {}
	for c in splittag:
		total_content[c[0]] = np.int_(c[1:])
	return total_content

## \From Anton

## Run through all reads and count all conversions

def errorRates(bamfilename,outfile,p_e_outfile):
	cellname='_'.join(bamfilename.split('/')[-1].split('_')[0:3])
	bamfile=pysam.AlignmentFile(bamfilename, 'rb')
	TC,SC,TC_ercc,SC_ercc=[],[],[],[]
	x,x_ercc=0,0
	for read in bamfile.fetch():
		if 'ERCC' in read.get_tag('XT'):
			x_ercc+=1
			TC_ercc.append(parseTCTag(read))
			SC_ercc.append(parseSCTag(read))
		else:
			x+=1
			TC.append(parseTCTag(read))
			SC.append(parseSCTag(read))
	tc=pd.DataFrame.from_dict(TC).sum(axis=0)
	sc=pd.DataFrame.from_dict(SC).sum(axis=0)
	scperc=sc.filter(regex='a')/tc['a']
	scperc=scperc.append(sc.filter(regex='c')/tc['c'])
	scperc=scperc.append(sc.filter(regex='g')/tc['g'])
	scperc=scperc.append(sc.filter(regex='t')/tc['t'])
	tc_ercc=pd.DataFrame.from_dict(TC_ercc).sum(axis=0)
	sc_ercc=pd.DataFrame.from_dict(SC_ercc).sum(axis=0)
	scperc_ercc=sc_ercc.filter(regex='a')/tc_ercc['a']
	scperc_ercc=scperc_ercc.append(sc_ercc.filter(regex='c')/tc_ercc['c'])
	scperc_ercc=scperc_ercc.append(sc_ercc.filter(regex='g')/tc_ercc['g'])
	scperc_ercc=scperc_ercc.append(sc_ercc.filter(regex='t')/tc_ercc['t'])
	dataArray={'ConversionType':scperc.keys(),'ConversionRate':scperc,'ConversionRate_ERCC':scperc_ercc,'Cell':[cellname]*scperc.shape[0],'ReadCount':[x]*scperc.shape[0],'ReadCount_ercc':[x_ercc]*scperc.shape[0]}
	df=pd.DataFrame(data=dataArray)
	if os.path.exists(outfile):
		print("Error Rate files are being overwritten...")
	df.to_csv(outfile)
	text_file=open(p_e_outfile,'w')
	text_file.write(str((scperc['cT']+scperc['gA'])/2))
	text_file.close()

if __name__ == "__main__":
	import os
	import pysam
	import pandas as pd
	import sys
	bamfilename = sys.argv[1]
	p_e_outfile = sys.argv[3]
	outfile = sys.argv[2]
	print('Calculating error rates for {}'.format(bamfilename))
	errorRates(bamfilename,outfile,p_e_outfile)

	#index outputed bamfile out of pipeline for lazyiness reason