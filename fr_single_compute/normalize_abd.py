import sys
import getopt
import os
import numpy as np
import pandas as pd


'''
    This is alpha_diversity.
    options:
    --abdf   <str> input file path
    --sname    <str> sample name on header
    --odir    <str> output directory
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

idir = '.'
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'sname=', 'odir='])
for op, arg in ops:
    if op == '--abdf':
        abd_path = arg
    if op == '--sname':
        sname = arg
    if op == '--odir':
        odir = arg



if not os.path.exists(odir):
    os.makedirs(odir)

abd = pd.read_csv(abd_path, sep='\t', header=None, index_col=0)
level = 's'
take_sp_list = []
rename_dict = {}
for idx in abd.index:
    last_l = idx.split('|')[-1]
    if last_l[0] == level:
        take_sp_list.append(idx)
        rename_dict[idx] = last_l
abd = abd.loc[take_sp_list, ]
abd.rename(index=rename_dict, inplace=True)
abd = abd.div(abd.sum(axis=0))
abd.rename(columns={1: sname}, inplace=True)

abd.to_csv(os.path.join(idir, '{}_relative_abundance.tsv'.format(sname)), sep='\t', header=True, index=True)


    
