from util import *
import sys
import subprocess
import getopt
import os

'''
    This is alpha_diversity.
    options:
    --abdf   <str> input file path
    --dsf    <str> input file from subplatform
    --ann    <str> input file of group info
    --groupid  <str> column name used for grouping, default: phenotype
    --method <str> method for alpha diversity, default: shannon. [shannon/simpson/invsimpson/ACE/Chao1/observedSpecies]
    --reads  <int> number of reads used to calculate richness, need with methods: ACE/Chao1/observedSpecies, default: 500000
    --testing <str> method for testing, default: t.test. [wilcox.test/t.test/kruskal.test/aov]
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)
script_path = os.path.join(python_dir, 'alpha.R')

ifile1 = ''
ifile2 = ''

ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'reads=', 'groupid=',"testing="])
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
    if op == '--reads':
        reads = arg
    if op == '--testing':
        testing = arg

merged_df = get_merged(ifile1,ifile2)
group = metadata2gf(metadata,groupid)
if not check_valid(group, merged_df):
    exit(2)
group_file = 'merged_input.group_info.tsv'
group.to_csv(group_file, sep='\t')


merged_df.to_csv('merged_input.abundance.all.tsv', sep='\t')
split_tax = tax_split(merged_df)
tax_list = [i for i in split_tax.keys()]

output1 = 'output.alpha_diversity.tsv'
output2 = 'output.pvalue.tsv'

for i in range(len(tax_list)):
    tax = tax_list[i]
    split_tax[tax].index.name = tax
    tax_abd = 'merged_input.abundance.' + tax +'.tsv'
    split_tax[tax].to_csv(tax_abd, sep='\t')
    if i == 0:
        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-r',reads,'-e',testing,'-a','F'],stdout=subprocess.PIPE)
    else:
        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-r',reads,'-e',testing,'-a','T'],stdout=subprocess.PIPE)
    if result.returncode > 0:
        print (result.stderr)
        exit(1)
