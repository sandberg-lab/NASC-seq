
#Static Baseline error rate at the moment
#p_e = 10e-4

##Pseudocode for sections 3.3,3.4 and 3.5 in GRAND-SLAM paper
##to be used together with a separate estimate for p_e as in 3.2

import pysam
import pandas as pd
import numpy as np
import scipy as sc
import pystan
from scipy.stats import binom
from scipy.special import beta
from scipy.stats import beta as beta_dist
import pickle
import stan_utility

def parseTCTag(read):
    tag = read.get_tag('TC')
    splittag = tag.split(';')
    total_content = {}
    for c in splittag:
        total_content[c[0]] = np.int_(c[1:])
    return total_content

def parseSCTag(read):
    tag = read.get_tag('SC')
    splittag = tag.split(';')
    specific_conversions = {}
    for c in splittag:
        specific_conversions[(c[0], c[1])] = np.int_(c[2:])
    return specific_conversions
def createAkng(Bamfile):
    A = {} #needs to be shrunk.. perhaps call g on demand
    for read in Bamfile:
        if read.get_tag('ST') == '+':
            k=parseSCTag(read)[('t','C')]
            n=parseTCTag(read)['t']
            g=read.get_tag('XT')

            if g in A:
                A[g][k,n] = A[g][k,n] + 1
            else:
                A[g] = np.zeros((300,300))
                A[g][k,n] = 1
        else:
            k=parseSCTag(read)[('a','G')]
            n=parseTCTag(read)['a']
            g=read.get_tag('XT')
            if g in A:
                A[g][k,n] = A[g][k,n] + 1
            else:
                A[g] = np.zeros((300,300))
                A[g][k,n] = 1
    return A

def createAkn(A):
    Akn = np.zeros(A[list(A.keys())[1]][1].shape) #does is remove the item?
    for g, Akng in A.items():
        Akn = Akn + Akng
    return Akn


def createMkn(Akn,p_e): #Left out from Akn
    M=np.zeros(Akn.shape)
    for n in range(Akn.shape[1]):
        for k in range(Akn.shape[0]):
            Ekn= np.sum(Akn[(k+1):,n])*binom.pmf(k,n,p_e)
            if Ekn > 0.01*Akn[k,n]:
                M[k,n]=1
    return M

def EstepAkn(Akn,Mkn,p_c): # Alters Akn in place - modify initial step
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            if Mkn[k,n]==1:
                num=0
                denom=0
                for kp in range(Mkn.shape[0]):
                    if Mkn[kp,n]==1:
                        num = num + binom.pmf(k,n,p_c)*Akn[kp,n]
                        denom = denom + binom.pmf(kp,n,p_c)
                Akn[k,n]=num/denom
    return Akn

def MstepP_c(Akn):
    num=0
    denom=0
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            num=num+k*Akn[k,n]
            denom=denom+n*Akn[k,n]
    p_c=num/denom
    return p_c

#Bisection-search for p_c  #Move code above into method to get local namespace
def estimateP_c(p_e,Bamfile):
    l=0
    r=1
    p_c0=(l+r)/2

    Akng=createAkng(Bamfile)
    Akn0=createAkn(Akng)

    p_c=p_c0
    Mkn=createMkn(Akn0,p_e)
    Akn=Akn0

    while r-l >= 10e-8:
        #print(r-l)
        Akn=EstepAkn(Akn,Mkn,p_c)
        p_c_old=p_c
        p_c=MstepP_c(Akn)
        if p_c < p_c_old:
            r=p_c
        else:
            l=p_c
    return p_c
## Code for determining Posterior pi_g (Unimodal Betadensities are the only cases possible by method in paper)
def createCounts(Bamfile):
    'Returns mode and parameters from Posterior Betaproportion of nascent reads'
    SC={}
    TC={}
    #pi_g={}

    #def logP(y,n,p_c,p_e):
    #    'logP (7) in paper, curried over intergration interval'
    #    def logPx(pg):
    #        return np.log((1-pg)*binom.pmf(y,n,p_e)+pg*binom.pmf(y,n,p_c))
    #    return logPx


    #def addPointwise(f,g):
    #    def summed(*args, **kwargs):
    #        return f(*args, **kwargs) + g(*args, **kwargs)
    #    return summed

    print('Adding counts for each read')
    for read in Bamfile.fetch():
        g = read.get_tag('XT')
        #print(g)
        if  read.get_tag('ST') == '+':
            y=parseSCTag(read)[('t','C')]
            n=parseTCTag(read)['t']
        else:
            y=parseSCTag(read)[('a','G')]
            n=parseTCTag(read)['a']

        #log-likelihood created as decorated functions per gene, adding a term for each read
        if g in SC:
            SC[g] = np.append(SC[g],y)
            TC[g] = np.append(TC[g],n)
        else:
            SC[g] = np.array([y])
            TC[g] = np.array([n])
    return SC,TC

def formatData(SC_g, TC_g):
    return dict(reads = len(SC_g), content = np.int_(TC_g), conversions = np.int_(SC_g), p_c = p_c, p_e = p_e)

def init_fun(dictionary):
    return [dictionary for i in np.arange(4)]

if __name__ == "__main__":
    import sys
    ifile = sys.argv[1]
    outdir = sys.argv[2]
    p_e = float(sys.argv[3])
    gene_list_file = sys.argv[4]
    stan_model = sys.argv[5]
    with open(gene_list_file) as f:
        gene_list = f.readlines()
    gene_list = [gene.strip() for gene in gene_list]
    Bam = pysam.AlignmentFile(ifile)#'P1_test_stranded.bam'
    model = stan_utility.compile_model(stan_model)
    try:
        print('Estimating p_c in {}'.format(ifile))
        p_c = estimateP_c(p_e,Bam)
        print('Estimated p_c is {}'.format(p_c))
    except RuntimeError:
        p_c = -1
        print('{} failed'.format(ifile))
    Bam = pysam.AlignmentFile(ifile)
    SC,TC = createCounts(Bam)
    params = {}
    meta_data = {}

    for gene in gene_list:
        if gene in SC.keys():
            print('Estimating {} in cell {}'.format(gene, ifile))
            data = formatData(SC[gene], TC[gene])
            meta_data[gene] = data
            fit =  model.sampling(data=data, seed=194838, init=init_fun(dict(alpha_logged=0, beta_logged=0, pi_g=0.5)))
            s = fit.summary()
            summary = pd.DataFrame(s['summary'], columns=s['summary_colnames'], index=s['summary_rownames'])
            params[gene] = summary
    params['meta_data'] = pd.DataFrame(meta_data)
    print('Wrote result of {} to {}'.format(ifile, outdir))
    pickle.dump(params, open(outdir, "wb"))
