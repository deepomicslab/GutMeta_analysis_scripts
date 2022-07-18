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
'''

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'diff.R')

ifile1 = ''
ifile2 = ''

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"adjust="])
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
    if op == '--adjust':
        adjust = arg


merged_df = get_merged(ifile1,ifile2)
group = metadata2gf(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = 'merged_input.group_info.tsv'
group.to_csv(group_file, sep='\t', na_rep='NA')



merged_df.to_csv('merged_input.abundance.all.tsv', sep='\t')
split_tax = tax_split(merged_df)
tax_list = [i for i in split_tax.keys()]

mean_cutoff = '0'
occ_cutoff = '0'
pvalue_cutoff = '0.05'

for i in range(len(tax_list)):
    tax = tax_list[i]
    split_tax[tax].index.name = tax
    new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
    split_tax[tax].index = new_tax
    tax_abd = 'merged_input.abundance.' + tax +'.tsv'
    split_tax[tax].to_csv(tax_abd, sep='\t')

    output = 'output.diff_testing.'+ method + '.' + tax 
    result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output,'-t',tax,'-m',method,'-a',adjust,'-p',pvalue_cutoff,'-e',mean_cutoff,'-c',occ_cutoff],stdout=subprocess.PIPE)
    if result.returncode > 0:
        print (result.stderr)
        exit(1)
    subprocess.run(["rm", tax_abd])
    
os.system("head -1 output.diff_testing.wilcox.test.s.stat.tsv|sed 's/^s/ID/' > output.taxonomy_tree.testing.tsv")
os.system("for i in `ls output.diff_testing.wilcox.test.*stat.tsv`; do sed 1d $i >> output.taxonomy_tree.testing.tsv;done")
os.system("rm output.diff_testing*")
