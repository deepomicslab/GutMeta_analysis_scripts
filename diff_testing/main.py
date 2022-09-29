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
    --method <str> method for testing [wilcox.test/t.test/kruskal.test/aov]
    --adjust <str> method for adjusting p-values, default: fdr. [fdr/BH/hochberg/holm/bonferroni/hommel]
    --mean_cutoff <float> threshold of mean value to filter taxonomy, at least one group should have a mean relative abundance greater than this threshold, defualt: 0
    --occ_cutoff <float> threshold of occurence to filter taxonomy, the taxonomy should at least exist in the threshold ratio of samples of certain group, defualt: 0.1
    --pvalue_cutoff <float> threshold of adjusted p-value as significant, defualt: 0.05
'''

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'diff.R')

ifile1 = ''
ifile2 = ''

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"adjust=","pvalue_cutoff=","occ_cutoff=","mean_cutoff="])
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
    if op == '--pvalue_cutoff':
        pvalue_cutoff = arg
    if op == '--occ_cutoff':
        occ_cutoff = arg
    if op == '--mean_cutoff':
        mean_cutoff = arg


merged_df = get_merged(ifile1,ifile2)
group = metadata2gf(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = 'merged_input.group_info.tsv'
group.to_csv(group_file, sep='\t', na_rep='NA')



merged_df.to_csv('merged_input.abundance.all.tsv', sep='\t')
split_tax = tax_split(merged_df)
tax_list = [i for i in split_tax.keys()]


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
        print ('ERROR:',result.stderr.decode())
        exit(1)
