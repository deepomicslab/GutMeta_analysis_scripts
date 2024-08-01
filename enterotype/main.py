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
    --method <str> method for visualization [PCA/PCoA]
    --maxk <int> maxinum number of clusters to find the optimal k. default:10
    --outdir <str> output directory, default: current directory
'''

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'enterotype.R')

ifile1 = ''
ifile2 = ''
outdir = '.'

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"maxk=", "outdir="])
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
    if op == '--maxk':
        max_k = arg
    if op == '--outdir':
        outdir = arg

merged_df = get_merged(ifile1,ifile2)
group = metadata2gf(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = os.path.join(outdir, 'merged_input.group_info.tsv')
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

    used_max_k = max_k
    if len(split_tax[tax].index) < int(max_k):
        used_max_k = len(split_tax[tax].index) - 1 
    if i != 0 :
        output = os.path.join(outdir, 'output.enterotype.' + fullname_dict[tax])
        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-o',output,'-t',tax,'-m',method,'-k',str(used_max_k)],stdout=subprocess.PIPE)
        print(result.stdout)
        if result.returncode > 0:
            print (result.stderr)
            exit(1)
    

for t in fullname_dict.values():
    if t == 'kingdom':
        continue
    coord_df = pd.read_csv(os.path.join(outdir, 'output.enterotype.{}.{}.coordinates.tsv'.format(t, method)), sep='\t', index_col=0, header=0)
    cluster_df = pd.read_csv(os.path.join(outdir, 'output.enterotype.' + t + '.cluster.tsv'), sep='\t', index_col=0, header=0)
    os.rename(os.path.join(outdir, 'output.enterotype.{}.{}.coordinates.tsv'.format(t, method)), os.path.join(outdir, 'output.enterotype.{}.coordinates.tsv'.format(t)))
    cluster_df['sample'] = cluster_df.index
    cluster_df['group'] = group.loc[cluster_df.index, groupid]
    cluster_df['enterotype'] = ['C{}'.format(x) for x in list(cluster_df['fac'])]
    cluster_df[['sample', 'group', 'enterotype']].to_csv(os.path.join(outdir, 'output.enterotype.{}.cluster.tsv'.format(t)), sep='\t', index=False)
    rename_dict  = {}
    for c in coord_df.columns:
        rename_dict[c] = c.replace('PC', 'Axis')
    coord_df = coord_df.rename(columns=rename_dict)
    coord_df.to_csv(os.path.join(outdir, 'output.enterotype.{}.coordinates.tsv'.format(t)), sep='\t', index=False)