from util import *
import sys
import subprocess
import getopt
import os
'''
    options:
    --abdf   <str> input file path
    --dsf    <str> input file from subplatform
    --ann    <str> input file of group info
    --groupid  <str> column name used for grouping, default: phenotype
'''
ifile1 = ''
ifile2 = ''

tax_level = 'genus'
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"tax_level="])
for op, arg in ops:
    if op == '--abdf':
        ifile1 = arg
    if op == '--dsf':
        ifile2 = arg
    if op == '--ann':
        metadata = arg
    if op == '--groupid':
        groupid = arg
    if op == '--tax_level':
        tax_level = arg



merged_df = get_merged(ifile1,ifile2)
group = filter_metadata_att(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
if 'Group' in group.columns:
    if groupid != 'Group':
        print('Warning: Attribute Group have already exist in metadata, no longer use choose groupid {}.'.format(groupid))
else:
    group = group.rename(columns={groupid:'Group'})
group_file = 'merged_input.metadata.tsv'
group.to_csv(group_file, sep='\t', na_rep='NA')


split_tax = tax_split(merged_df)

tax_dict = {'species': 's', 'genus': 'g', 'class': 'c', 'order': 'o', 'phylum': 'p'}
tax = tax_dict[tax_level]

split_tax[tax].index.name = tax
new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
split_tax[tax].index = new_tax

df=split_tax[tax].T
s = df.sum()
tmp = df[s.sort_values(ascending=False).index[:20]].T
tmp.loc['Other'] = 100 - tmp.sum()
tax_abd = 'merged_input.abundance.top20.' + tax +'.tsv'
tmp.to_csv(tax_abd, sep='\t')
