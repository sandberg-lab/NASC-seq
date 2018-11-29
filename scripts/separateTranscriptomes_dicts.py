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
        print(key) 
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
def separateTranscriptomes(pkl, total_df):
    no_expr = (total_df == 0)

    ratio_nascent = getMode(pkl)
    ratio_nascent = ratio_nascent.reindex(total_df.index, axis=0)
    ratio_nascent = ratio_nascent.reindex(total_df.columns, axis=1)
    ratio_nonnascent = 1-ratio_nascent

    nascent_df = ratio_nascent.multiply(total_df, axis='index')
    nonnascent_df = ratio_nonnascent.multiply(total_df, axis='index')

    nascent_df = nascent_df.fillna(-1)
    nonnascent_df = nonnascent_df.fillna(-1)

    nascent_df = nascent_df.add(no_expr, axis='index')
    nonnascent_df = nonnascent_df.add(no_expr, axis='index')

    nascent_df = nascent_df.replace(-1, np.nan)
    nonnascent_df = nonnascent_df.replace(-1, np.nan)
    return nascent_df, nonnascent_df
if __name__ == '__main__':
    pkl_list_file = sys.argv[1]
    counttable = sys.argv[2]
    outname = sys.argv[3]
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
    print('Loading {}'.format(counttable))
    total_df = pd.read_csv(counttable, index_col=0)
    nascent_df, nonnascent_df = separateTranscriptomes(dict_list, total_df)
    nascenttablename = '{}_{}.csv'.format(outname,'nascentTable')
    print('Saving nascent count table to {}'.format(nascenttablename))
    nascent_df.to_csv(nascenttablename)
    nonnascenttablename = '{}_{}.csv'.format(outname,'noNnascentTable')
    print('Saving non-nascent count table to {}'.format(nonnascenttablename))
    nonnascent_df.to_csv(nonnascenttablename)