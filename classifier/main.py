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
    --pvalue_cutoff <str> threshold of adjusted p-value as significant, defualt: 0.05
    --occ_cutoff <str> threshold of occurence to filter taxonomy, the taxonomy should at least exist in the threshold ratio of samples of certain group, defualt: 0.1
    --mean_cutoff <str> threshold of mean value to filter taxonomy, at least one group should have a mean relative abundance greater than this threshold, defualt: 0
    --tax_level <str> taxonomy level to find biomarker, default: species. [species/genus]
    --scale <boolean> scale abundance with log or not, defualt: T [T:True/F:False]
    --cv_n <int> CV fold. default = 10 [2-10]
    --rep_n <int> replicates number for cross-validation. default = 5 [2-10]
    --outdir <str> output directory, default: current directory
'''

script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'diff.R')
script_path2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'random_forest.R')

ifile1 = ''
ifile2 = ''
outdir = '.'

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'groupid=',"adjust=","pvalue_cutoff=","occ_cutoff=","mean_cutoff=","tax_level=","scale=","cv_n=","rep_n=", 'outdir='])
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
    if op == '--tax_level':
        tax_level = arg
    if op == '--scale':
        scale = arg
        if scale == 'T':
            scale = 'True'
        elif scale == 'F':
            scale = 'False'
    if op == '--cv_n':
        cv_n = arg
    if op == '--rep_n':
        rep_n = arg
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

tax_dict = {'species': 's', 'genus': 'g'}

tax = tax_dict[tax_level]
split_tax[tax].index.name = tax
new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
split_tax[tax].index = new_tax
tax_abd = os.path.join(outdir, 'merged_input.abundance.' + tax +'.tsv')
split_tax[tax].to_csv(tax_abd, sep='\t')

output = os.path.join(outdir, 'output.diff_testing.'+ method + '.' + tax )
result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output,'-t',tax,'-m',method,'-a',adjust,'-p',pvalue_cutoff,'-e',mean_cutoff,'-c',occ_cutoff],stdout=subprocess.PIPE)
if result.returncode > 0:
    print (result.stderr)
    exit(1)

pass_abd = os.path.join(outdir, output + '.abd.pass.tsv')
rf_output = os.path.join(outdir, 'output.random_forest.' + method + '.' + tax_level + '.pass.' + rep_n + '_' + cv_n)
result = subprocess.run(['Rscript',script_path2,'-i',pass_abd,'-g',group_file,'-o',rf_output,'-s',scale,'-r',rep_n,'-c',cv_n],stdout=subprocess.PIPE)

    
if result.returncode > 0:
    print (result.stderr)
    exit(1)
