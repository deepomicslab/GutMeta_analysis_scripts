import sys
import subprocess
import getopt
import os
import pandas as pd
import copy
from scipy.stats import fisher_exact

def overlap(range1, range2):
    if range1[0] > range2[1] or range1[1] < range2[0]:
        return False
    else:
        return True

def search_event(scaffold, db, range):
    scaffold_df = db[db['Chr'] == scaffold]
    valid_idx = []
    for idx in scaffold_df.index:
        if overlap(range, [scaffold_df.loc[idx, 'Start'], scaffold_df.loc[idx, 'End']]):
            valid_idx.append(idx)
    valid_df = scaffold_df.loc[valid_idx, ]
    return valid_df
    

def search_row(idx, df, db_dir, fr_size):
    tmp = '{}.VFDB.tsv'
    row = df.loc[idx, ]
    # for recipient
    recipient = row['recipient']
    chrom = recipient.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    db = pd.read_csv(ifile, header=0, index_col=None, sep='\t')
    recipient_range = [max(0, row['insert_locus']-fr_size), row['insert_locus']+fr_size]
    recipient_df = search_event(recipient, db, recipient_range)
    # for donor
    donor = row['donor']
    range = [max(0, row['delete_start'] - fr_size), row['delete_end'] + fr_size]
    chrom = donor.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    db = pd.read_csv(ifile, header=0, index_col=None, sep='\t')
    donor_df = search_event(donor, db, range)
    return recipient_df, donor_df 

def enrichment(metadata, result_anno, groupid, outdir):
    pheno_set = list(set(metadata[groupid]))
    g1 = pheno_set[0]
    g2 = pheno_set[1]
    vf_df = result_anno[(result_anno['recipient_VF_n']>0) | (result_anno['donor_VF_n']>0)]
    cate_df = pd.DataFrame(columns=[g1, g2])
    for idx in vf_df.index:
        recipient_VF_category = vf_df.loc[idx, 'recipient_VF_category'].split(';')
        donor_VF_category = vf_df.loc[idx, 'donor_VF_category'].split(';')
        all_cates = recipient_VF_category+donor_VF_category
        for cate in all_cates:
            if cate == 'NA':
                continue
            if cate not in cate_df.index:
                cate_df.loc[cate, g1] = 0
                cate_df.loc[cate, g2] = 0
            cate_df.loc[cate, metadata.loc[vf_df.loc[idx, 'sample'], groupid]] += 1
    pvalue_reformat = pd.DataFrame(columns=['group1', 'group2', 'category', 'g1_in_category', 'g1_total', 'g2_in_category', 'g2_total', 'pvalue', 'odds_ratio'])
    g1_total = cate_df[g1].sum()
    g2_total = cate_df[g2].sum()
    for cate in cate_df.index:
        a = cate_df.loc[cate, g1]
        b = cate_df.loc[cate, g2]
        c = g1_total - a
        d = g2_total - b
        if a+b == 0 or c+d == 0 or a+c == 0 or b+d == 0:
            pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, 'NA', 'NA']
            continue
        oddsratio, pvalue = fisher_exact([[a, b], [c, d]])
        pvalue_reformat.loc[cate, ] = [g1, g2, cate, a, g1_total, b, g2_total, pvalue, oddsratio]
    return pvalue_reformat

'''
    options:
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --fr_size    <int> flanking region size
    --ann   <str> metadata file
    --groupid  <str> group id
    --outdir    <str> output directory
'''
ops, args = getopt.getopt(sys.argv[1:], '', ['db_dir=', 'hgt=', 'fr_size=', 'ann=', 'groupid=', 'outdir='])
db_file = '/data2/platform/gutmeta_v2_platform/Database/genome/DB.genome_annotation'
fr_size = 1000
outdir = '.'
for op, arg in ops:
    if op == '--db_dir':
        db_idir = arg
    if op == '--hgt':
        infile = arg
    if op == '--fr_size':
        fr_size = int(arg)
    if op == '--ann':
        infile2 = arg
    if op == '--groupid':
        groupid = arg
    if op == '--outdir':
        outdir = arg

if not os.path.exists(outdir):
    os.makedirs(outdir)
#infile = '../../HGT_demo_file/SAMEA3449210.event_output.csv'
#db_file = '../../HGT_demo_file/HGT/DB.HGT_clusters.annotated.tsv'


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

MGE_result = pd.DataFrame()
result_anno = pd.DataFrame(columns=['id', 'sample', 'recipient_VF_n', 'recipient_VF_category', 'recipient_VF_list', 'donor_VF_n', 'donor_VF_category', 'donor_VF_list', 'recipient', 'insert_locus', 'donor', 'delete_start', 'delete_end', 'reverse_flag'])
for idx in df.index:
    recipient_df, donor_df = search_row(idx, df, db_idir, fr_size)
    id = 'HGT_c{}'.format(idx+1)
    sample = df.loc[idx, 'sample']
    recipient_MGE_n = recipient_df.shape[0]
    recipient_MGE_category = ';'.join(recipient_df['Category'])
    recipient_MGE_list = ';'.join(recipient_df['Name'])
    if recipient_MGE_n == 0:
        recipient_MGE_list = 'NA'
        recipient_MGE_category = 'NA'
    donor_MGE_n = donor_df.shape[0]
    donor_MGE_category = ';'.join(recipient_df['Category'])
    donor_MGE_list = ';'.join(recipient_df['Name'])
    if donor_MGE_n == 0:
        donor_MGE_list = 'NA'
        donor_MGE_category = 'NA'
    recipient = df.loc[idx, 'recipient']
    insert_locus = df.loc[idx, 'insert_locus']
    donor = df.loc[idx, 'donor']
    delete_start = df.loc[idx, 'delete_start']
    delete_end = df.loc[idx, 'delete_end']
    reverse_flag = df.loc[idx, 'reverse_flag']
    result_anno.loc[len(result_anno), ] = [id, sample, recipient_MGE_n, recipient_MGE_category, recipient_MGE_list, donor_MGE_n, donor_MGE_category, donor_MGE_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]
    #result_anno.iloc[len(result_anno), ] = [id, sample, recipient_HGTC_n, recipient_HGTC_list, donor_HGTC_n, donor_HGTC_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]
    merge_df = pd.concat([recipient_df, donor_df], ignore_index=True)
    if len(MGE_result)==0:
        MGE_result = copy.deepcopy(merge_df)
    else:
        MGE_result = pd.concat([MGE_result, merge_df], ignore_index=True)
MGE_result.drop_duplicates(inplace=True)
MGE_result.to_csv(os.path.join(outdir, 'output.VF_comparison.VF.tsv'), index=False, sep='\t')
result_anno.to_csv(os.path.join(outdir, 'output.VF_comparison.annotated.tsv'), index=False, sep='\t')
pvalue_df = enrichment(metadata, result_anno, groupid, outdir)
pvalue_df.to_csv(os.path.join(outdir, 'output.VF_comparison.pvalue.tsv'), index=False, sep='\t')