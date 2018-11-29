
import pandas as pd
import pickle
import sys
import os
import numpy as np
def iter_over_dicts(dict_list):
    for dictionary in dict_list:
        for key, value in dictionary.items():
            yield key, value
def getParam_dict(dict_list,par, estim):
    df = pd.DataFrame()
    for key, pkl in iter_over_dicts(dict_list):
        value_list = []
        for gene in pkl.keys():
            value_list.append(pkl[gene][estim][par])
        series = pd.Series(value_list, index=pkl.keys())
        series.name = key
        df = df.append(series)
    return df
def getMode(dict_list):
    alpha = getParam_dict(dict_list, 'alpha', 'mean')
    beta = getParam_dict(dict_list, 'beta', 'mean')
    pi_g_mode = (alpha - 1)/(alpha + beta - 2)
    pi_g_mode = pi_g_mode.where(~(alpha < 1).rmul(beta >= 1), other=0)
    pi_g_mode = pi_g_mode.where(~(alpha  >= 1).rmul(beta < 1), other=1)
    return pi_g_mode
if __name__ == '__main__':
    pkl_list_file = sys.argv[1]
    outname = sys.argv[2]
    print('Loading pickles')
    with open(pkl_list_file) as f:
        pickle_list = f.readlines()
        print(pickle_list)
    pickle_list = [pickle.strip() for pickle in pickle_list]
    print(pickle_list)
    dict_list = []
    for pkl in pickle_list:
        print('Loading in {}'.format(pkl))
        dict = pickle.load(open(pkl, 'rb'))
        dict_list.append(dict)
    df = getMode(dict_list)
    tablename = '{}_{}.csv'.format(outname,'mode')
    print('Saving mode to {}'.format(tablename))
    df.to_csv(tablename)
