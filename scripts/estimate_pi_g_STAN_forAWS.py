
#Static Baseline error rate at the moment
#p_e = 10e-4

##Pseudocode for sections 3.3,3.4 and 3.5 in GRAND-SLAM paper
##to be used together with a separate estimate for p_e as in 3.2


import pandas as pd
import numpy as np
import logging
logger = logging.getLogger("pystan")
# add root logger (logger Level always Warning)
# not needed if PyStan already imported
logger.addHandler(logging.NullHandler())
logger_path = "/dev/null"
fh = logging.FileHandler(logger_path, encoding="utf-8")
fh.setLevel(logging.INFO)
import pystan
import pickle
from joblib import Parallel, delayed

def formatData(SC_g, TC_g, p_c_cell, p_e_cell):
    if len(SC_g) > 10000:
        random_reads = np.random.choice(len(SC_g), 1000, replace=False)
        SC_g = SC_g[random_reads]
        TC_g = TC_g[random_reads]
    return {'reads' : len(SC_g), 'content' : np.int_(TC_g), 'conversions' : np.int_(SC_g), 'p_c' : p_c_cell, 'p_e' : p_e_cell}
    
def init_fun(dictionary):
    return [dictionary for i in np.arange(4)]

def estim_pi_g(gene_key,gene, model):
            print('Estimating {}'.format(gene_key))
            data = formatData(gene['SC'], gene['TC'], gene['p_c'], gene['p_e'])
            if data['reads'] >= 16:
                fit =  model.sampling(data=data, n_jobs=1, seed=194838, init=init_fun({'alpha_logged':0, 'beta_logged':0, 'pi_g':0.5}))
                s = fit.summary()
                summary = pd.DataFrame(s['summary'], columns=s['summary_colnames'], index=s['summary_rownames'])
                return summary, gene_key
            else:
                return None
def iter_over_dicts(dict_list):
    for dictionary in dict_list:
        for key, value in dictionary.items():
            yield key, value

if __name__ == "__main__":
    import sys
    pickle_list_file = sys.argv[1]
    outdir = sys.argv[2]
    stanFile = sys.argv[3]
    model = pystan.StanModel(file=stanFile,verbose=True)
    with open(pickle_list_file) as f:
        pickle_list = f.readlines()
        print (pickle_list)
    pickle_list = [pikl.strip() for pikl in pickle_list]
    print (pickle_list)
    dict_list = []
    for pkl in pickle_list:
        print('Loading in {}'.format(pkl))
        dic = pickle.load(open(pkl, 'rb'))
        dict_list.append(dic)

    params = Parallel(n_jobs=72, verbose = 3)(delayed(estim_pi_g)(gene_key, gene, model) for (gene_key,gene) in iter_over_dicts(dict_list))
    master_dict = {}
    for df, key in params:
        if df is not None:
            both_keys = key.split('-')
            if len(both_keys) == 2:
                gene_key = both_keys[0]
                cell_key = both_keys[1]
            else:
                gene_key = ''.join(both_keys[:-1])
                cell_key = both_keys[-1]
            if cell_key in master_dict:
                master_dict[cell_key][gene_key] = df
            else:
                master_dict[cell_key] = {}
                master_dict[cell_key][gene_key] = df

    print('Wrote result to'.format(outdir))
    pickle.dump(master_dict, open(outdir, "wb"))