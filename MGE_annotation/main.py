import sys
import subprocess
import getopt
import os
import pandas as pd

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
    tmp = '{}.ImmeDB.tsv'
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


'''
    options:
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --fr_size    <int> flanking region size
    --outdir    <str> output directory
'''
ops, args = getopt.getopt(sys.argv[1:], '', ['db_dir=', 'hgt=', 'fr_size=', 'outdir='])
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
    if op == '--outdir':
        outdir = arg


#infile = '../../HGT_demo_file/SAMEA3449210.event_output.csv'
#db_file = '../../HGT_demo_file/HGT/DB.HGT_clusters.annotated.tsv'


df = pd.read_csv(infile, header=0, index_col=None)
df.rename(columns={'receptor':'recipient'}, inplace=True)

result_anno = pd.DataFrame(columns=['id', 'sample', 'recipient_MGE_n', 'recipient_MGE_category', 'recipient_MGE_list', 'donor_MGE_n', 'donor_MGE_category', 'donor_MGE_list', 'recipient', 'insert_locus', 'donor', 'delete_start', 'delete_end', 'reverse_flag'])
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
merge_df.to_csv(os.path.join(outdir, 'output.MGE_annotation.MGE.tsv'), index=False, sep='\t')
result_anno.to_csv(os.path.join(outdir, 'output.MGE_annotation.annotated.tsv'), index=False, sep='\t')
