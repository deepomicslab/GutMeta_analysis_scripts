import sys
from util import *

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

taxonomy_fullname = 'kingdom/phylum/class/order/family/genus/species'
fullname_dict = {}
for t in taxonomy_fullname.split('/'):
    fullname_dict[t[0]] = t

diversity = pd.read_csv('output.alpha_diversity.tsv', sep='\t', index_col=0, header=0)
pvalue = pd.read_csv('output.pvalue.tsv', sep='\t', index_col=0, header=0)
reformat_diversity = pd.DataFrame(columns=['sample', 'group', 'taxonomy_level', 'alpha_diversity'])
phenos = pd.read_csv('merged_input.group_info.tsv', sep='\t', index_col=0, header=0)
for s in diversity.columns:
    pheno = phenos.loc[s, groupid]
    for taxa_level in diversity.index:
        idx = len(reformat_diversity)
        reformat_diversity.loc[idx, 'sample'] = s
        reformat_diversity.loc[idx, 'taxonomy_level'] = fullname_dict[taxa_level]
        reformat_diversity.loc[idx, 'alpha_diversity'] = diversity.loc[taxa_level, s]
        reformat_diversity.loc[idx, 'group'] = pheno
reformat_diversity.to_csv('output.alpha_diversity.results.tsv', sep='\t', index=False)

pheno_set = list(set(phenos[groupid]))
g1 = pheno_set[0]
g2 = pheno_set[1]
pvalue_reformat = pd.DataFrame(index=list(fullname_dict.values()), columns=['group1', 'group2', 'g1_mean', 'g1_variance', 'g2_variance', 'g2_mean', 'g1/g2', 'enrich', 'pvalue'])
for taxa_level in pvalue.index:
    pvalue_reformat.loc[fullname_dict[taxa_level], 'group1'] = g1
    pvalue_reformat.loc[fullname_dict[taxa_level], 'group2'] = g2
    pvalue_reformat.loc[fullname_dict[taxa_level], 'g1_mean'] = reformat_diversity[(reformat_diversity['group'] == g1)^(reformat_diversity['taxonomy_level'] == fullname_dict[taxa_level])]['alpha_diversity'].mean()
    pvalue_reformat.loc[fullname_dict[taxa_level], 'g1_variance'] = reformat_diversity[(reformat_diversity['group'] == g1)^(reformat_diversity['taxonomy_level'] == fullname_dict[taxa_level])]['alpha_diversity'].var()
    pvalue_reformat.loc[fullname_dict[taxa_level], 'g2_mean'] = reformat_diversity[(reformat_diversity['group'] == g2)^(reformat_diversity['taxonomy_level'] == fullname_dict[taxa_level])]['alpha_diversity'].mean()
    pvalue_reformat.loc[fullname_dict[taxa_level], 'g2_variance'] = reformat_diversity[(reformat_diversity['group'] == g2)^(reformat_diversity['taxonomy_level'] == fullname_dict[taxa_level])]['alpha_diversity'].var()
    pvalue_reformat.loc[fullname_dict[taxa_level], 'g1/g2'] = pvalue_reformat.loc[fullname_dict[taxa_level], 'g1_mean'] / pvalue_reformat.loc[fullname_dict[taxa_level], 'g2_mean']
    pvalue_reformat.loc[fullname_dict[taxa_level], 'pvalue'] = pvalue.loc[taxa_level, 'pvalue']
    if (pvalue_reformat.loc[fullname_dict[taxa_level], 'g1_mean'] > pvalue_reformat.loc[fullname_dict[taxa_level], 'g2_mean']):
        pvalue_reformat.loc[fullname_dict[taxa_level], 'enrich'] = g1
    else:
        pvalue_reformat.loc[fullname_dict[taxa_level], 'enrich'] = g2

pvalue_reformat.to_csv('output.alpha_diversity.pvalue.tsv', sep='\t')