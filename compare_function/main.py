import sys
import subprocess
import getopt
import os
import pandas as pd
import gzip
import COG_enrichment as ce
import kegg_enrichment as ke
import pickle
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection as fdr


def gz2df(ifile):
    with gzip.open(ifile, 'rb') as f:
        df = pd.read_csv(f, header=None, index_col=None, sep='\t')
    return df

def extract_line(content):
    ko_list = []
    cog_list = []
    cog_cate_list = []
    # check kegg
    if 'KEGG=' in content:
        kegg = content.split('KEGG=')[1].strip().split(';')[0]
        kos = kegg.split(',')
        #print('kos', kos)
        for ko in kos:
            if ko.startswith('ko:'):
                ko = ko.split('ko:')[-1]
                ko_list.append(ko)
    
    if 'db_xref=' in content:
        cog = content.split('db_xref=')[1].strip().split(';')[0]
        cogs = cog.split(',')
        for cog in cogs:
            if cog.startswith('COG:'):
                cog = cog.split('COG:')[-1]
                cog_list.append(cog)
    if 'COG=' in content:
        cog_cate = content.split('COG=')[1].strip().split(';')[0]
        cog_cate = cog_cate.split(',')[0]
        for cc in cog_cate:
            if cc != '-':
                cog_cate_list += [char for char in cc]
    
    return ko_list, cog_list, cog_cate_list

def extract(df):
    ko_list = []
    cog_list = []
    cog_cate_list = []
    for idx in df.index:
        kos, cogs, cog_cates = extract_line(df.loc[idx, 8])
        ko_list += kos
        cog_list += cogs
        cog_cate_list += cog_cates

    return ko_list, cog_list, cog_cate_list

def overlap(range1, range2):
    if range1[0] > range2[1] or range1[1] < range2[0]:
        return False
    else:
        return True

# styp = recipient or donor
def search_event(scaffold, db, range):
    scaffold_df = db[db[0] == scaffold]
    valid_idx = []
    for idx in scaffold_df.index:
        if overlap(range, [scaffold_df.loc[idx, 3], scaffold_df.loc[idx, 4]]):
            valid_idx.append(idx)
    valid_df = scaffold_df.loc[valid_idx, ]
    return valid_df
    
def search_row(idx, df, db_dir, fr_size):
    tmp = '{}.gff.gz'
    row = df.loc[idx, ]
    # for recipient
    recipient = row['recipient']
    chrom = recipient.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    db = gz2df(ifile)
    recipient_range = [max(0, row['insert_locus']-fr_size), row['insert_locus']+fr_size]
    recipient_df = search_event(recipient, db, recipient_range)
    
    # for donor
    donor = row['donor']
    range = [max(0, row['delete_start'] - fr_size), row['delete_end'] + fr_size]
    chrom = donor.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    db = gz2df(ifile)
    donor_df = search_event(donor, db, range)
    return recipient_df, donor_df

def get_backgroud(scaffold, db, range_excluded):
    scaffold_df = db[db[0] == scaffold]
    valid_idx = []
    for idx in scaffold_df.index:
        if not overlap(range_excluded, [scaffold_df.loc[idx, 3], scaffold_df.loc[idx, 4]]):
            valid_idx.append(idx)
    valid_df = scaffold_df.loc[valid_idx, ]
    return valid_df


'''
    options:
    --db_dir   <str> input dir of KEGG database
    --hgt    <str> input file of HGT output
    --fr_size    <int> flanking region size
    --ko_pathway_dict <str> input file of ko_pathway_dict
    --ann   <str> metadata file
    --groupid  <str> group id
    --table <str> input file of pathway table
    --outdir <str> output dir
'''

ops, args = getopt.getopt(sys.argv[1:], '', ['db_dir=', 'hgt=', 'fr_size=', 'ko_pathway_dict=', 'outdir=', 'ann=', 'groupid=', 'table='])
db_file = '/data2/platform/gutmeta_v2_platform/Database/genome/DB.genome_annotation'
ko_pathway_dict = '/data2/platform/gutmeta_v2_platform/Database/function_db/ko_pathway_dict.pickle'
fr_size = 1000
outdir = '.'
for op, arg in ops:
    if op == '--db_dir':
        db_idir = arg
    if op == '--hgt':
        infile = arg
    if op == '--fr_size':
        fr_size = int(arg)
    if op == '--ko_pathway_dict':
        pfile = arg
    if op == '--outdir':
        outdir = arg
    if op == '--ann':
        infile2 = arg
    if op == '--table':
        name_table = arg
    if op == '--groupid':
        groupid = arg

if not os.path.exists(outdir):
    os.makedirs(outdir)

pathway_table = pd.read_csv(name_table, sep='\t', header=0, index_col=0)
with open(pfile, 'rb') as f: 
    ko_pathway_dict = pickle.load(f)
df = pd.read_csv(infile, header=0, index_col=None)
df.rename(columns={'receptor':'recipient'}, inplace=True)
metadata = pd.read_csv(infile2, header=0, index_col=0, sep='\t')
if len(metadata[groupid].unique()) != 2:
    print('Error: the column {} does not have exact 2 level.'.format(groupid))
    exit(1)
hgt_slist = list(set(df['sample']))
g_slist = list(metadata.index)
valid = True
for s in list(hgt_slist):
    if s not in g_slist:
        print('Error: group information of sample {} in HGT event dose NOT exist.'.format(s))
        valid = False
if not valid:
    exit(2)

result_anno = pd.DataFrame(columns=['id', 'sample', 'phenotype', 'recipient_KEGG_n', 'recipient_KEGG_list', 'recipient_COG_n', 'recipient_COG_list',
                                    'donor_KEGG_n', 'donor_KEGG_list', 'donor_COG_n', 'donor_COG_list',
                                    'recipient', 'insert_locus', 'donor', 'delete_start', 'delete_end', 'reverse_flag'])
cog_relation = {}
related_set_dict = {}
for idx in df.index:
    recipient_df, donor_df = search_row(idx, df, db_idir, fr_size)
    id = 'HGT_c{}'.format(idx+1)
    sample = df.loc[idx, 'sample']
    phenotype = metadata.loc[sample, groupid]
    ko_list, cog_list, cog_cate = extract(recipient_df)
    recipient_KEGG_n = len(set(ko_list))
    recipient_COG_n = len(set(cog_list))
    recipient_KEGG_list = ';'.join(set(ko_list))
    recipient_COG_list = ';'.join(set(cog_list))
    if recipient_KEGG_n == 0:
        recipient_KEGG_list = 'NA'
    if recipient_COG_n == 0:
        recipient_COG_list = 'NA'
    if phenotype not in related_set_dict.keys():
        related_set_dict[phenotype] = {'ko':set(), 'cog':set(), 'cog_cate':[]}
    related_set_dict[phenotype]['ko'] = related_set_dict[phenotype]['ko'].union(set(ko_list))
    related_set_dict[phenotype]['cog'] = related_set_dict[phenotype]['cog'].union(set(cog_list))
    related_set_dict[phenotype]['cog_cate'] += related_set_dict[phenotype]['cog_cate']
    for i, cog in enumerate(cog_list):
        cog_relation[cog] = cog_cate[i]
    ko_list, cog_list, cog_cate = extract(donor_df)
    donor_KEGG_n = len(set(ko_list))
    donor_COG_n = len(set(cog_list))
    donor_KEGG_list = ';'.join(set(ko_list))
    donor_COG_list = ';'.join(set(cog_list))
    if donor_KEGG_n == 0:
        donor_KEGG_list = 'NA'
    if donor_COG_n == 0:
        donor_COG_list = 'NA'
    related_set_dict[phenotype]['ko'] = related_set_dict[phenotype]['ko'].union(set(ko_list))
    related_set_dict[phenotype]['cog'] = related_set_dict[phenotype]['cog'].union(set(cog_list))
    related_set_dict[phenotype]['cog_cate'] += related_set_dict[phenotype]['cog_cate']
    for i, cog in enumerate(cog_list):
        cog_relation[cog] = cog_cate[i]
    recipient = df.loc[idx, 'recipient']
    insert_locus = df.loc[idx, 'insert_locus']
    donor = df.loc[idx, 'donor']
    delete_start = df.loc[idx, 'delete_start']
    delete_end = df.loc[idx, 'delete_end']
    reverse_flag = df.loc[idx, 'reverse_flag']
    result_anno.loc[len(result_anno), ] = [id, sample, phenotype, recipient_KEGG_n, recipient_KEGG_list, recipient_COG_n, recipient_COG_list, donor_KEGG_n, donor_KEGG_list, donor_COG_n, donor_COG_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]
result_anno.to_csv(os.path.join(outdir, 'output.functional_annotation.annotated.tsv'), index=False, sep='\t')

# kegg to pathway count
cate_df = pd.DataFrame(columns=related_set_dict.keys())
for pheno in related_set_dict.keys():
    for ko in related_set_dict[pheno]['ko']:
        if ko in ko_pathway_dict.keys():
            pathways = ko_pathway_dict[ko]
            for pathway in pathways:
                if not pathway.startswith('map'):
                    continue
                if pathway not in cate_df.index:
                    cate_df.loc[pathway, pheno] = 0
                    cate_df.fillna(0, inplace=True)
                cate_df.loc[pathway, pheno] += 1
cate_df.fillna(0, inplace=True)

pheno_set = list(set(metadata[groupid]))
g1 = pheno_set[0]
g2 = pheno_set[1]
pvalue_reformat = pd.DataFrame(columns=['group1', 'group2', 'category', 'g1_in_category', 'g1_total', 'g2_in_category', 'g2_total', 'pvalue', 'odds_ratio'])
g1_total = len(related_set_dict[g1]['ko'])
g2_total = len(related_set_dict[g2]['ko'])
valid_cate = []
for cate in cate_df.index:
    a = cate_df.loc[cate, g1]
    b = cate_df.loc[cate, g2]
    c = g1_total - a
    d = g2_total - b
    if a+b==0 or c+d==0 or a+c==0 or b+d==0:
        pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, 'NA', 'NA']
        continue
    oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
    if b*c == 0:
        oddsratio = 'NA'
    pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, pvalue, oddsratio]
    
    valid_cate.append(cate)
padj = fdr(pvalue_reformat.loc[valid_cate, 'pvalue'].tolist(), 0.05)[1]
for i, cate in enumerate(valid_cate):
    pvalue_reformat.loc[cate, 'padj'] = padj[i]
    
for idx in pvalue_reformat.index:
    pname, fc, sc = ke.get_pathway_name_class_static(idx, pathway_table)
    pvalue_reformat.loc[idx, 'pathway_name'] = pname
    pvalue_reformat.loc[idx, 'first_class'] = fc
    pvalue_reformat.loc[idx, 'second_class'] = sc
pvalue_reformat.sort_values(by=['pvalue'], ascending=True).to_csv(os.path.join(outdir, 'output.function_comparison.pvalue.KEGG.tsv'), index=False, sep='\t')

# cog to pathway count
cate_df = pd.DataFrame(columns=related_set_dict.keys())
for pheno in related_set_dict.keys():
    for cog  in related_set_dict[pheno]['cog']:
        cogc = cog_relation[cog]
        if cogc not in cate_df.index:
            cate_df.loc[cogc, pheno] = 0
            cate_df.fillna(0, inplace=True)
        cate_df.loc[cogc, pheno] += 1
cate_df.fillna(0, inplace=True)

pheno_set = list(set(metadata[groupid]))
g1 = pheno_set[0]
g2 = pheno_set[1]
pvalue_reformat = pd.DataFrame(columns=['group1', 'group2', 'category', 'g1_in_category', 'g1_total', 'g2_in_category', 'g2_total', 'pvalue', 'odds_ratio'])
g1_total = len(related_set_dict[g1]['cog'])
g2_total = len(related_set_dict[g2]['cog'])
valid_cate = []
for cate in cate_df.index:
    a = cate_df.loc[cate, g1]
    b = cate_df.loc[cate, g2]
    c = g1_total - a
    d = g2_total - b
    if a+b==0 or c+d==0 or a+c==0 or b+d==0:
        pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, 'NA', 'NA']
        continue
    oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
    pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, pvalue, oddsratio]
    valid_cate.append(cate)
padj = fdr(pvalue_reformat.loc[valid_cate, 'pvalue'].tolist(), 0.05)[1]
for i, cate in enumerate(valid_cate):
    pvalue_reformat.loc[cate, 'padj'] = padj[i]

COG_dict, COG_profile_dict = ce.get_COG_dict()
for idx in pvalue_reformat.index:
    pvalue_reformat.loc[idx, 'category'] = COG_dict[idx]
    pvalue_reformat.loc[idx, 'profile'] = COG_profile_dict[idx]
pvalue_reformat.sort_values(by=['pvalue'], ascending=True).to_csv(os.path.join(outdir, 'output.function_comparison.pvalue.COG.tsv'), index=False, sep='\t')
