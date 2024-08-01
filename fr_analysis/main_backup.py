import sys
import getopt
import os
import numpy as np
import pandas as pd

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
    This is alpha_diversity.
    options:
    --abdf   <str> input file path of related abundance
    --gcn_d    <str> input GCN
    --top_n   <int> top n species / <float> top n ratio # 0 or negative or not input means all
    --odir    <str> output directory
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

odir = '.'
top_n = 0
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'gcn_d=', 'top_n=', 'odir='])
for op, arg in ops:
    if op == '--abdf':
        abd_path = arg
    if op == '--gcn_d':
        d_path = arg
    if op == '--top_n':
        top_n = float(arg)
    if op == '--odir':
        odir = arg



if not os.path.exists(odir):
    os.makedirs(odir)

abd = pd.read_csv(abd_path, sep='\t', header=0, index_col=0)
d_df = pd.read_csv(d_path, sep='\t', index_col=0, header=0)
nfr_df = pd.DataFrame(index=abd.columns, columns=['nFR'])
sname  = abd.columns[0]
nFR, fr_df, profile = nfr(d_df, abd, sname)
nfr_df.loc[sname, 'nFR'] = nFR
row_index, col_index = np.tril_indices(len(fr_df), k=0)
fr_df.values[row_index, col_index] = 0
edge_df = fr_df.stack().reset_index()
edge_df = edge_df[edge_df[0] != 0]
edge_df.sort_values(by=0, ascending=False, inplace=True)
if top_n >= 1:
    edge_df = edge_df.head(int(top_n))
    top_n = '{:.0f}'.format(top_n)
elif top_n > 0:
    edge_df = edge_df.head(int(top_n*len(edge_df)))
else:
    top_n = 'all'

with open(os.path.join(odir, '{}.nFR.txt'.format(sname)), 'w') as f:
    f.write(str(nFR))
edge_df.to_csv(os.path.join(odir, '{}.nFR.top[{}].edge.tsv'.format(sname, top_n)), sep='\t', header=False, index=False)


    
