import pysam
import numpy as np
import pandas as pd
import pickle
from scipy.stats import binom

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
def createAkng(Bamfile, cell_id):
    A = {} #needs to be shrunk.. perhaps call g on demand
    SC = {}
    TC = {}
    for read in Bamfile.fetch():
        g='{}_{}'.format(read.get_tag('XT'), cell_id)
        if read.get_tag('ST') == '+':
            k=parseSCTag(read)[('t','C')]
            n=parseTCTag(read)['t']

            if g in A:
                SC[g] = np.append(SC[g],k)
                TC[g] = np.append(TC[g],n)
                A[g][k,n] = A[g][k,n] + 1
            else:
                SC[g] = np.array([k])
                TC[g] = np.array([n])
                A[g] = np.zeros((300,300))
                A[g][k,n] = 1
        else:
            k=parseSCTag(read)[('a','G')]
            n=parseTCTag(read)['a']

            if g in A:
                SC[g] = np.append(SC[g],k)
                TC[g] = np.append(TC[g],n)
                A[g][k,n] = A[g][k,n] + 1
            else:
                SC[g] = np.array([k])
                TC[g] = np.array([n])
                A[g] = np.zeros((300,300))
                A[g][k,n] = 1
    return A, SC, TC

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
def estimateP_c(p_e,Bamfile, cell_id):
    l=0
    r=1
    p_c0=(l+r)/2

    Akng, SC, TC=createAkng(Bamfile, cell_id)
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
    return SC, TC, p_c, p_e
if __name__ == "__main__":
    import sys
    ifile = sys.argv[1]
    outdir = sys.argv[2]
    p_e = float(sys.argv[3])
    cell_id = sys.argv[4]
    p_c = -1
    if len(sys.argv) == 6:
        p_c = float(sys.argv[5])
    Bam = pysam.AlignmentFile(ifile)#'P1_test_stranded.bam'
    if p_c == -1:

        try:
            print('Estimating p_c in {}'.format(ifile))
            SC, TC, p_c, p_e = estimateP_c(p_e,Bam, cell_id)
            print('Estimated p_c is {}'.format(p_c))
            concatenated_dict = {}
            for g in SC.keys():
                concatenated_dict[g] = {'SC': SC[g], 'TC': TC[g], 'p_c': p_c, 'p_e':p_e}
            print('Saving prepared data to {}'.format(outdir))
            pickle.dump(concatenated_dict, open(outdir, 'wb'))
        except RuntimeError:
            p_c = -1
            print('{} failed'.format(ifile))
    else:
        concatenated_dict = {}
        Akng, SC, TC=createAkng(Bam, cell_id)
        for g in SC.keys():
            concatenated_dict[g] = {'SC': SC[g], 'TC': TC[g], 'p_c': p_c, 'p_e':p_e}
        print('Saving prepared data to {}'.format(outdir))
        pickle.dump(concatenated_dict, open(outdir, 'wb'))
