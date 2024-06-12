import os
import pandas as pd
import getopt
import sys

def rename_x(df):
    rename_dict = {}
    for idx in df.index:
        rename_dict[idx] = idx.replace('x__', 'k__')
    df = df.rename(index=rename_dict)
    return df

def merge_dir(idir):
    merged_abd = pd.DataFrame()
    for file in os.listdir(idir):
        df = pd.read_csv(os.path.join(idir, file), sep='\t', header=None, index_col=0)
        sname = file.split('.')[0]
        df.dropna(how='all', axis=0)
        df.dropna(how='all', axis=1)
        df = df.loc[:, (df != 0).any(axis=0)]
        df = df.loc[(df != 0).any(axis=1), :]
        col_sum = df.sum()
        normalized_df = df.div(col_sum)
        normalized_df = normalized_df*7
        normalized_df = normalized_df.rename(columns={1: sname})
        abd = normalized_df[sname]
        merged_abd = pd.concat([merged_abd, abd], axis=1)
    merged_abd.fillna(0, inplace=True)
    return merged_abd

def merge_slist(idir, slist):
    merged_abd = pd.DataFrame()
    for file in os.listdir(idir):
        sname = file.split('.')[0]
        if sname not in slist:
            continue
        df = pd.read_csv(os.path.join(idir, file), sep='\t', header=None, index_col=0)
        df.dropna(how='all', axis=0)
        df.dropna(how='all', axis=1)
        df = df.loc[:, (df != 0).any(axis=0)]
        df = df.loc[(df != 0).any(axis=1), :]
        col_sum = df.sum()
        normalized_df = df.div(col_sum)
        normalized_df = normalized_df*7
        normalized_df = normalized_df.rename(columns={1: sname})
        abd = normalized_df[sname]
        merged_abd = pd.concat([merged_abd, abd], axis=1)
    merged_abd.fillna(0, inplace=True)
    return merged_abd

def merge_slistfile(idir, slistfile):
    with open(slistfile) as f:
        slist = f.read().split('\n')
    return merge_slist(idir, slist)

def merge_abd(idir, slist=None):
    if slist is None:
        return merge_dir(idir)
    if os.path.isfile(slist):
        return merge_slistfile(idir, slist)
    return merge_slist(idir, slist)

'''
    This is used to merge bracken result.
    options:
    --idir: input directory <str>
    --slist: sample list file  <str> (optional)
    --opath: output path <str>
'''

ops, args = getopt.getopt(sys.argv[1:], '', ['idir=', 'slist=', 'opath='])
slist = None
for op, arg in ops:
    if op == '--idir':
        idir = arg
    if op == '--slist':
        slist = arg
    if op == '--opath':
        opath = arg

merged_abd = merge_abd(idir, slist)
merge_abd = rename_x(merged_abd)
merged_abd.to_csv(opath, sep='\t')