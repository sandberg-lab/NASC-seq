
#Static Baseline error rate at the moment
#p_e = 10e-4

##Pseudocode for sections 3.3,3.4 and 3.5 in GRAND-SLAM paper
##to be used together with a separate estimate for p_e as in 3.2

import pysam
import pandas as pd
import numpy as np
import scipy as sc
from scipy.stats import binom
from scipy.special import beta
from scipy.stats import beta as beta_dist
import warnings
warnings.filterwarnings("ignore")

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
    for read in Bamfile.fetch():
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
    Akn = np.zeros((300,300))
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
def createLL_g(Bamfile,p_c,p_e):
    'Returns mode and parameters from Posterior Betaproportion of nascent reads'
    P_c={}
    P_e={}
    alpha_dict={}
    beta_dict={}
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
        if g in P_c:
            P_c[g] = np.append(P_c[g],binom.pmf(y,n,p_c))
            P_e[g] = np.append(P_e[g],binom.pmf(y,n,p_e))
        else:
            P_c[g] = np.array([binom.pmf(y,n,p_c)])
            P_e[g] = np.array([binom.pmf(y,n,p_e)])
            alpha_dict[g] = 1  #prior has 1
            beta_dict[g]=1
    return P_c,P_e,alpha_dict,beta_dict

def calcPg(P_cg,P_eg,alpha_g,beta_g, g):
        'Returns mode and parameters from Posterior Betaproportion of nascent reads'
        from scipy.special import beta
        from scipy.stats import beta as beta_fun
        def logLL( pc_array, pe_array):
            def LogLP(pg):
                return np.sum(np.log(pg*pc_array + (1-pg)*pe_array))
            return LogLP
        def BisectionForAdaptiveIntInterval(f,mode): #f to be integrated
            l1=0
            r1=mode
            l2=mode
            r2=1
            l=l1
            h=r2
            i = 0
            while r1-l1 > 10e-8 and r2-l2 >10e-8:
                if f(l) < 10e-3*f(mode):
                    l=(l1+r1)/2
                    l1=l
                if f(l) > 10e-3*f(mode):
                    l=(l1+r1)/2
                    r1=l
                if f(h) < 10e-3*f(mode):
                    h=(l2+r2)/2
                    r2=h
                if f(h) > 10e-3*f(mode):
                    h=(l2+r2)/2
                    l2=h
                if i > 1000:
                    return np.nan,np.nan
                i += 1
            return l,h
        def betadist(x,a,b):#(pg,alpha_g,beta_g):
            'density of standard beta distribution'
            return (x**(a-1)*(1-x)**(b-1))/beta(a,b)

        def betadistc(a,b):#(pg,alpha_g,beta_g):
            'parameters of density of standard beta distribution'
            def density(x):
                return (x**(a-1)*(1-x)**(b-1))/beta(a,b)
            return density
            #beta is the beta function int_{0,1}[x**(a-1)*(1-x)**(b-1)]
        pg0=0.5 # initial guess for mode
        k=20 #n.o. intervals in numerial integration
        def fixPointwise(f,g):
            def summed(*args, **kwargs):
                return (f(*args, **kwargs) - g(*args, **kwargs))**2
            return summed
        P = {}
        P[g] = logLL(P_cg, P_eg)
        #change to map evaluation over genes g
        #def dens(g,P):
        #    'Currying on gene'
        #    def density(x):
        #        pg = x[0]
        #        alpha_g = x[1]
        #        beta_g = x[2]
        #        return (P[g](pg)+np.log(beta_fun.pdf(pg,alpha_g,beta_g)))
        #    return density
        def dens(g,P):
            'Currying on gene'
            def density(x):
                pg = x[0]
                alpha_g = x[1]
                beta_g = x[2]
                return -(P[g](pg)+np.log(beta_fun.pdf(pg,alpha_g,beta_g)))
            return density
        x0=[pg0,alpha_g,beta_g]
        fun = dens(g,P)

        result=sc.optimize.minimize(fun,x0,bounds=((0,1),(0,1e5),(0,1e5)),method='L-BFGS-B', options={'maxiter':100})
        max_val=result.fun
        maxpg=result.x[0]
        alpha_g=result.x[1]
        beta_g=result.x[2]
        success = result.success
        #change to map evaluation over genes g
        def denspg(g,P,alpha_g,beta_g):
            'Currying on optimized values per gene'
            def density(pg):
                return np.exp(P[g](pg)+np.log(beta_fun.pdf(pg,alpha_g,beta_g)))
            return density
        def denspg1(g,P,alpha_g,beta_g):
            'Currying on optimized values per gene'
            def density(pg):
                return -np.exp(P[g](pg)+np.log(beta_fun.pdf(pg,alpha_g,beta_g)))
            return density

        l,h = BisectionForAdaptiveIntInterval(denspg(g,P,alpha_g,beta_g),maxpg)
        if np.isnan(l):
            return g, np.nan,np.nan,np.nan

        # Compute integral by Trapezoidrule and optimize to find best betadistribution fit (actually a bit shady- we force the posterior to be a Beta)
        # need to curry betadist
        def L2target(k,l,h,g,P,denspg):
            'Target for Bayesian optimization'
            def trapetsRule(k,l,h,denspg):
                A=(denspg(l)+denspg(h))/2
                for i in range(k-2):
                    A+= denspg(l+(h-l)/k*(i+1))
                return (h-l)/k*A
              #closed form not likely as important as they claim in paper
            def valuefunc(x):
                alpha_g=x[0]
                beta_g=x[1]
                return trapetsRule(k,l,h,fixPointwise(denspg(g,P,alpha_g,beta_g),betadistc(alpha_g,beta_g)))
            return valuefunc

        #Curry dens-betadistr into trapezoidrule
        x0=[alpha_g,beta_g]
        fun = L2target(k,l,h,g,P,denspg)
        result=sc.optimize.minimize(fun,x0,bounds=((0,1e5),(0,1e5)),method='L-BFGS-B', options={'maxiter':100})
        #result.fun=max_val
        alpha_g=result.x[0]
        beta_g=result.x[1]
        return g,alpha_g,beta_g, success

def ComputeOpt(P_c,P_e,alpha_dict,beta_dict):
    params = [calcPg(P_c[g], P_e[g], alpha_dict[g], beta_dict[g], g) for g in P_c.keys()]
    return params


if __name__ == "__main__":
    import sys
    ifile = sys.argv[1]
    outdir = sys.argv[2]
    p_e = float(sys.argv[3])
    Bam = pysam.AlignmentFile(ifile)#'P1_test_stranded.bam'
    try:
        print('Estimating p_c in {}'.format(ifile))
        p_c = estimateP_c(p_e,Bam)
        print(p_c)
    except RuntimeError:
        p_c = -1
        print('{} failed'.format(ifile))
    Bam = pysam.AlignmentFile(ifile)
    print('Adding logP for each read for {}'.format(ifile))
    P_c,P_e,alpha_dict,beta_dict = createLL_g(Bam,p_c,p_e)
    params = ComputeOpt(P_c,P_e,alpha_dict,beta_dict)
    pg_df = pd.DataFrame(params, columns=['gene', 'alpha', 'beta', 'optim_success']).dropna()
    pg_df = pg_df.set_index('gene')
    mode_dict = {}
    for gene, par in pg_df.iterrows():
        if  par['alpha'] < 1 and par['beta'] < 1:
            mode_dict[gene] = np.nan
        elif np.absolute(1-par['alpha']) < 0.1 and np.absolute(1-par['beta']) < 0.1:
            mode_dict[gene] = np.nan
        elif par['alpha'] > 1 and par['beta'] > 1:
            mode_dict[gene] = (par['alpha'] -1)/(par['alpha']+ par['beta'] -2)
        elif par['alpha'] < 1 and par['beta'] > 1:
            mode_dict[gene] = 0
        elif par['alpha'] > 1 and par['beta'] < 1:
            mode_dict[gene] = 1
    mode_series = pd.Series(mode_dict, index=mode_dict.keys())
    mode_series.name = 'pi_g'
    pg_df = pg_df.join(mode_series)
    print('saving to {}'.format(outdir))
    pg_df.to_csv(outdir)
