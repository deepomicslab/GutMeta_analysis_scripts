import sys
import getopt
import os
import numpy as np
import pandas as pd
from util import *
from scipy.stats import ttest_ind, kruskal, f_oneway
from scipy.stats import mannwhitneyu

def nfr(d_df, profile, sname):
    sp_list = list(set(list(profile.index)).intersection(set(list(d_df.index))))
    n = len(sp_list)
    corr = np.ones(shape=(n, n)) - d_df.loc[sp_list, sp_list].values
    np.fill_diagonal(corr, 0)
    a = np.array(profile.loc[sp_list, sname])
    inter_matrix = np.dot(a.reshape(len(a), 1),a.reshape(1, len(a)))
    np.fill_diagonal(inter_matrix, 0)
    td = np.sum(inter_matrix)/2
    fr = np.sum(np.multiply(inter_matrix, corr))/2
    fr_df = pd.DataFrame(np.multiply(inter_matrix, corr), index=sp_list, columns=sp_list)
    profile = profile.loc[sp_list, sname]
    if td == 0:
        return 0
    return fr/td, fr_df, profile


'''
    This is nfr comparison.
    options:
    --abdf   <str> input file path of related abundance
    --gcn_d    <str> input GCN distance
    --ann    <str> input file of group info
    --groupid  <str> column name used for grouping, default: phenotype
    --method <str> method for testing [wilcox.test/t.test/kruskal.test/aov]
    --odir    <str> output directory
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

odir = '.'
ifile1 = ''
ifile2 = ''
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'gcn_d=', 'ann=', 'odir=', 'groupid=', 'method='])
for op, arg in ops:
    if op == '--abdf':
        ifile1 = arg
    if op == '--gcn_d':
        d_path = arg
    if op == '--ann':
        metadata = arg
    if op == '--odir':
        odir = arg
    if op == '--groupid':
        groupid = arg
    if op == '--method':
        method = arg

if not os.path.exists(odir):
    os.makedirs(odir)


merged_df = get_merged(ifile1,ifile2)
group = metadata2gf(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = os.path.join(odir, 'merged_input.group_info.tsv')
group.to_csv(group_file, sep='\t', na_rep='NA')
split_tax = tax_split(merged_df)
sp_df = split_tax['s']

rename_dict = {}
for idx in sp_df.index:
    last_l = idx.split('|')[-1]
    rename_dict[idx] = last_l
sp_df.rename(index=rename_dict, inplace=True)
sp_df = sp_df.div(sp_df.sum(axis=0))

d_df = pd.read_csv(d_path, sep='\t', index_col=0, header=0)
nfr_df = pd.DataFrame(index=sp_df.columns, columns=['group', 'nFR'])

for sname in sp_df.columns:
    nFR, fr_df, profile = nfr(d_df, sp_df, sname)
    nfr_df.loc[sname, 'nFR'] = nFR
    nfr_df.loc[sname, 'group'] = group.loc[sname, groupid]

group = pd.read_csv(group_file, sep='\t', index_col=0)
phenos = list(set(group[groupid]))
g1 = phenos[0]
g2 = phenos[1]

cols = [
    'group1',
    'group2',
    'g1_mean',
    'g2_mean',
    'g1_variance',
    'g2_variance',
    'g1_n',
    'g2_n',
    'p_value',
    'g1/g2', 
    'enriched']



p_df = pd.DataFrame(columns=cols)
p_df.loc[0, 'group1'] = g1
p_df.loc[0, 'group2'] = g2
g1_v = nfr_df.loc[group[group[groupid] == g1].index, 'nFR'].values.astype(float)
g2_v = nfr_df.loc[group[group[groupid] == g2].index, 'nFR'].values.astype(float)
p_df.loc[0, 'g1_mean'] = g1_v.mean()
p_df.loc[0, 'g2_mean'] = g2_v.mean()
if p_df.loc[0, 'g1_mean'] > p_df.loc[0, 'g2_mean']:
    p_df.loc[0, 'enriched'] = g1
else:
    p_df.loc[0, 'enriched'] = g2
p_df.loc[0, 'g1/g2'] = p_df.loc[0, 'g1_mean']/p_df.loc[0, 'g2_mean']
p_df.loc[0, 'g1_variance'] = g1_v.var()
p_df.loc[0, 'g2_variance'] = g2_v.var()
#p_df.loc[0, 'g1_occ'] = len(g1_v[g1_v > 0])/len(g1_v)
#p_df.loc[0, 'g2_occ'] = len(g2_v[g2_v > 0])/len(g2_v)
p_df.loc[0, 'g1_n'] = len(g1_v)
p_df.loc[0, 'g2_n'] = len(g2_v)

if method == 't.test':
    p_df.loc[0, 'p_value'] = ttest_ind(g1_v, g2_v)[1]
elif method == 'wilcox.test':
    p_df.loc[0, 'p_value'] = mannwhitneyu(g1_v, g2_v)[1]
elif method == 'kruskal.test':
    p_df.loc[0, 'p_value'] = kruskal(g1_v, g2_v)[1]
elif method == 'aov':
    p_df.loc[0, 'p_value'] = f_oneway(g1_v, g2_v)[1]
else:
    print('Error: method not supported.')
    exit(1)

p_df.to_csv(os.path.join(odir, 'output.FR_comparison.testing.tsv'), sep='\t', index=False)
nfr_df.to_csv(os.path.join(odir, 'output.FR_comparison.nFR.tsv'), sep='\t', index=True)