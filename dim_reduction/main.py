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
    --method <str> method for dimension reduction, default: PCA. [PCA/ICA/tSNE/UMAP/NMDS/PCoA]
    --dist <str> method for distance method, work for NMDS/PCoA default: bray. [bray/euclidean/jaccard/manhattan]
    --axis <int> number of axis in the output, not work for tSNE, default: 10.
    --outdir <str> output directory, default: current directory
'''

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reduce_dim.R')

ifile1 = ''
ifile2 = ''
outdir = '.'

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"dist=","axis=", 'outdir='])
for op, arg in ops:
    if op == '--abdf':
        ifile1 = arg
    if op == '--dsf':
        ifile2 = arg
    if op == '--ann':
        metadata = arg
    if op == '--groupid':
        groupid = arg
    if op == '--method':
        method = arg
    if op == '--dist':
        dist_method = arg
    if op == '--axis':
        axis_n = arg
    if op == '--outdir':
        outdir = arg

merged_df = get_merged(ifile1,ifile2)
group = filter_metadata_att(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = os.path.join(outdir, 'merged_input.metadata.tsv')
group.to_csv(group_file, sep='\t', na_rep='NA')


merged_df.to_csv(os.path.join(outdir, 'merged_input.abundance.all.tsv'), sep='\t')
split_tax = tax_split(merged_df)
tax_list = [i for i in split_tax.keys()]

taxonomy_fullname = 'kingdom/phylum/class/order/family/genus/species'
fullname_dict = {}
for t in taxonomy_fullname.split('/'):
    fullname_dict[t[0]] = t

for i in range(len(tax_list)):
    tax = tax_list[i]
    split_tax[tax].index.name = tax
    new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
    split_tax[tax].index = new_tax
    tax_abd = os.path.join(outdir, 'merged_input.abundance.' + tax +'.tsv')
    split_tax[tax].to_csv(tax_abd, sep='\t')


    if i != 0 and (method == 'tSNE' or len(split_tax[tax].index) > int(axis_n)):
        output = os.path.join(outdir, 'output.dimension_reduction.{}.{}.tsv'.format(method,fullname_dict[tax]))
        #output = 'output.' + method + '.reduce_dim.' + tax + '.tsv'
        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-o',output,'-t',tax,'-m',method,'-r',axis_n,'-d',dist_method],stdout=subprocess.PIPE)
        if result.returncode > 0:
            print (result.stderr)
            exit(1)
