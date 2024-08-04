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
        if scale == 'T' or scale == 'True':
            scale = 'True'
        elif scale == 'F'or scale == 'False':
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
rename_dict = {}
for idx in merged_df.index:
    rename_dict[idx] = idx.replace('-', '_')
merged_df.rename(index=rename_dict, inplace=True)
merged_df.to_csv(os.path.join(outdir, 'merged_input.abundance.all.tsv'), sep='\t')
split_tax = tax_split(merged_df)

tax_dict = {'species': 's', 'genus': 'g'}

tax = tax_dict[tax_level]
split_tax[tax].index.name = tax
new_tax = [(tax_id.split('|')[-1]) for tax_id in list(split_tax[tax].index)]
split_tax[tax].index = new_tax
tax_abd = os.path.join(outdir, 'merged_input.abundance.tsv')
split_tax[tax].to_csv(tax_abd, sep='\t')

output = os.path.join(outdir, 'output.diff_testing')
result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output,'-t',tax,'-m',method,'-a',adjust,'-p',pvalue_cutoff,'-e',mean_cutoff,'-c',occ_cutoff],stdout=subprocess.PIPE)
if result.returncode > 0:
    print (result.stderr)
    exit(1)

pass_abd = os.path.join(outdir, output + '.abd.pass.tsv')
rf_output = os.path.join(outdir, 'temp.classifier')
result = subprocess.run(['Rscript',script_path2,'-i',pass_abd,'-g',group_file,'-o',rf_output,'-s',scale,'-r',rep_n,'-c',cv_n],stdout=subprocess.PIPE)

    
if result.returncode > 0:
    print (result.stderr)
    exit(1)

group = pd.read_csv(group_file, sep='\t', index_col=0)
phenos = list(set(group[groupid]))
g1 = phenos[0]
g2 = phenos[1]

col_rename = {
    'mean({})'.format(g1): 'g1_mean',
    'mean({})'.format(g2): 'g2_mean',
    'sd({})'.format(g1): 'g1_variance',
    'sd({})'.format(g2): 'g2_variance',
    'occ-rate({})'.format(g1): 'g1_occ',
    'occ-rate({})'.format(g2): 'g2_occ',
    'occ-n({})'.format(g1): 'g1_n',
    'occ-n({})'.format(g2): 'g2_n',
    'pvalue': 'p_value',
    'qvalue': 'p_adj',
}
taxonomy_fullname = 'genus/species'
fullname_dict = {}
for t in taxonomy_fullname.split('/'):
    fullname_dict[t[0]] = t

p_df_new = pd.DataFrame(columns=(['taxonomy_level', 'group1', 'group2']+list(col_rename.values())+['g1/g2', 'enriched']))
p_pass_df_new = pd.DataFrame(columns=(['taxonomy_level', 'group1', 'group2']+list(col_rename.values())+['g1/g2', 'enriched']))
abd_df_pass_new = pd.DataFrame(columns=group.index)
tax = tax_level
output_p = os.path.join(outdir, 'output.diff_testing.stat.tsv')
output_pass_a = os.path.join(outdir, 'output.diff_testing.abd.pass.tsv')
abd_df_pass = pd.read_csv(output_pass_a, sep='\t', index_col=0)
p_df = pd.read_csv(output_p, sep='\t', index_col=0)
p_df.rename(columns=col_rename, inplace=True)
p_df['taxonomy_level'] = tax_dict[tax]
p_df['g1/g2'] = p_df['g1_mean']/p_df['g2_mean']
p_df['g1_variance'] = p_df['g1_variance']**2
p_df['g2_variance'] = p_df['g2_variance']**2
p_df['group1'] = g1
p_df['group2'] = g2
for t in p_df.index:
    p_df_new.loc[t, :] = p_df.loc[t, list(p_df_new.columns)]
    if t in abd_df_pass.index:
        p_pass_df_new.loc[t, :] = p_df_new.loc[t, :]
        abd_df_pass_new.loc[t, ] = abd_df_pass.loc[t, abd_df_pass.columns]
p_df_new.to_csv(os.path.join(outdir, 'output.classifier.testing.tsv'), sep='\t')
p_pass_df_new.to_csv(os.path.join(outdir, 'output.classifier.testing.pass.tsv'), sep='\t')
# abd_df_pass_new.to_csv(os.path.join(outdir, 'output.classifier.testing.pass.abundance.tsv'), sep='\t')


# rename 
os.rename(os.path.join(outdir, 'temp.classifier.cross_validation.marker.predict.in.train.tsv'), os.path.join(outdir, 'output.classifier.prediction.tsv'))
os.rename(os.path.join(outdir, 'temp.classifier.train.auc_info.txt'), os.path.join(outdir, 'temp.classifier.note.tsv'))

importance_df = pd.read_csv(os.path.join(outdir, 'temp.classifier.feature.importance.tsv'), index_col=0, sep='\t')
with open(os.path.join(outdir, 'temp.classifier.cross_validation_pick.txt')) as fp:
    pick_list = fp.read().strip().split(',')
for sp in pick_list:
    importance_df.loc[sp, 'Selected'] = 'True'
importance_df.fillna('False', inplace=True)
importance_df.to_csv(os.path.join(outdir, 'output.classifier.feature_importance.tsv'), sep='\t')