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
    --method <str> method for beta diversity calculation [bray/euclidean/jaccard/manhattan/weighted_unifrac/unweighted_unifrac]
    --reads  <int> number of reads used to calculate richness, need with methods: ACE/Chao1/observedSpecies, default: 500000
    --testing <str> method for testing, default: t.test. [wilcox.test/t.test/kruskal.test/aov]
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)
script_path = os.path.join(python_dir, 'beta.R')

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


for i in range(len(tax_list)):
    tax = tax_list[i]
    split_tax[tax].index.name = tax
    new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
    split_tax[tax].index = new_tax
    tax_abd = 'merged_input.abundance.' + tax +'.tsv'
    split_tax[tax].to_csv(tax_abd, sep='\t')

    output1 = 'output.' + method + '.beta_diversity.' + tax
    output2 = 'output.' + method + '.pvalue.' + tax + '.tsv'

    if (method == "weighted_unifrac" or method == "unweighted_unifrac"):
        tree = 'merged_input.newick_tree.' + tax + '.txt'
        tree_f = open(tree,'w')
        tree_f.write(map2newick(split_tax[tax], 0.5))
        tree_f.close()

        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-r',reads,'-e',testing,'-c',tree],stdout=subprocess.PIPE)
    else:
        result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-e',testing],stdout=subprocess.PIPE)
    if result.returncode > 0:
        print (result.stderr)
        exit(1)

phenos = pd.read_csv('merged_input.group_info.tsv', sep='\t', index_col=0, header=0)
pheno_set = list(set(phenos[groupid]))
g1 = pheno_set[0]
g2 = pheno_set[1]
taxonomy_fullname = 'kingdom/phylum/class/order/family/genus/species'
fullname_dict = {}
for t in taxonomy_fullname.split('/'):
    fullname_dict[t[0]] = t

matrix_df_new = pd.DataFrame(columns=['sample_pair', 'group', 'taxonomy_level', 'beta_diversity'])
for t in fullname_dict.keys():
    taxonomy_fullname = fullname_dict[t]
    infile = 'output.{}.beta_diversity.{}.matrix.tsv'.format(method, t)
    matrix_df = pd.read_csv(infile, sep='\t', index_col=0, header=0)
    nsample = matrix_df.shape[0]
    samples = matrix_df.index
    matrix_df = matrix_df.loc[samples, samples]
    sample1_list = []
    sample2_list = []
    values = []
    for i in range(nsample):
        sample1_list += [samples[i]] * (nsample - i - 1)
        sample2_list += list(samples[i+1:])
        values += list(matrix_df.iloc[i, i+1:])
    pheno_list1 = [phenos.loc[s, groupid] for s in sample1_list]
    pheno_list2 = [phenos.loc[s, groupid] for s in sample2_list]
    ppair_list = ['{}:{}'.format(pheno_list1[i], pheno_list2[i]) for i in range(len(pheno_list1))]
    spair_list = ['{}:{}'.format(sample1_list[i], sample2_list[i]) for i in range(len(sample1_list))]
    matrix_df_new_t = pd.DataFrame({'sample_pair': spair_list, 'group': ppair_list, 'taxonomy_level': taxonomy_fullname, 'beta_diversity': values})
    matrix_df_new = pd.concat([matrix_df_new, matrix_df_new_t])
matrix_df_new.to_csv('output.beta_diversity.results.tsv ', sep='\t', index=False)

p_df_new = pd.DataFrame(columns=['taxonomy_level', 'group1', 'group2', 'g1_n', 'g2_n', 'g1_occ', 'g2_occ', 'g1_mean', 'g1_variance', 'g2_variance', 'g2_mean', 'g1/g2', 'enrich', 'pvalue'])
nsample = phenos.shape[0]
nsample = (nsample - 1) * nsample / 2
for t in fullname_dict.keys():
    taxonomy_fullname = fullname_dict[t]
    infile = 'output.{}.pvalue.{}.tsv'.format(method, t)
    p_df = pd.read_csv(infile, sep='\t', index_col=0, header=0)
    for ppair in p_df.index:
        g1 = ppair.split(':')[0].replace('-', ':')
        g2 = ppair.split(':')[1].replace('-', ':')
        samples1 = matrix_df_new[(matrix_df_new['group'] == g1)&(matrix_df_new['taxonomy_level'] == taxonomy_fullname)]['beta_diversity']
        samples2 = matrix_df_new[(matrix_df_new['group'] == g2)&(matrix_df_new['taxonomy_level'] == taxonomy_fullname)]['beta_diversity']
        g1_mean = samples1.mean()
        g2_mean = samples2.mean()
        g1_variance = samples1.var()
        g2_variance = samples2.var()
        g1_g2 = g1_mean / g2_mean
        if g1_g2 > 1:
            enrich = g1
        else:
            enrich = g2
        pvalue = p_df.loc[ppair, 'p-value']
        idx = p_df_new.shape[0]
        p_df_new.loc[idx, ] = [taxonomy_fullname, g1, g2, len(samples1), len(samples2), len(samples1)/nsample, len(samples2)/nsample, g1_mean, g1_variance, g2_variance, g2_mean, g1_g2, enrich, pvalue]
p_df_new.to_csv('output.beta_diversity.pvalue.tsv', sep='\t', index=False)
            

