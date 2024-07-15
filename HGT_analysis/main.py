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

# styp = recipient or donor
def search_event(scaffold, db, range, stype):
    scaffold_df = db[db['{}_scaffold'.format(stype)] == scaffold]
    valid_idx = []
    if stype == 'donor':
        for idx in scaffold_df.index:
            if overlap(range, [scaffold_df.loc[idx, 'donor_bp1'], scaffold_df.loc[idx, 'donor_bp2']]):
                valid_idx.append(idx)
        valid_df = scaffold_df.loc[valid_idx, ]
    else:
        valid_df = scaffold_df[(scaffold_df['recipient_bp'] >= range[0]) & (scaffold_df['recipient_bp'] <= range[1])]
    return valid_df

def search_row(idx, df, db, fr_size):
    row = df.loc[idx, ]
    # for recipient
    recipient = row['recipient']
    recipient_range = [max(0, row['insert_locus']-fr_size), row['insert_locus']+fr_size]
    recipient_df = search_event(recipient, db, recipient_range, 'recipient')
    # for donor
    donor = row['donor']
    range = [max(0, row['delete_start'] - fr_size), row['delete_end'] + fr_size]
    donor_df = search_event(donor, db, range, 'donor')
    return recipient_df, donor_df 


'''
    options:
    --db   <str> input file path of HGT cluster database
    --hgt    <str> input file of HGT output
    --fr_size    <int> flanking region size
    --outdir    <str> output directory
'''
outdir = '.'
ops, args = getopt.getopt(sys.argv[1:], '', ['db=', 'hgt=', 'fr_size=', 'outdir='])
db_file = '/data2/platform/gutmeta_v2_platform/Database/HGT/DB.HGT_clusters.annotated.tsv'
fr_size = 1000
for op, arg in ops:
    if op == '--db':
        db_file = arg
    if op == '--hgt':
        infile = arg
    if op == '--fr_size':
        fr_size = int(arg)
    if op == '--outdir':
        outdir = arg


#infile = '../../HGT_demo_file/SAMEA3449210.event_output.csv'
#db_file = '../../HGT_demo_file/HGT/DB.HGT_clusters.annotated.tsv'

df = pd.read_csv(infile, header=0, index_col=None)
db = pd.read_csv(db_file, header=0, index_col=0, sep='\t')
df.rename(columns={'receptor':'recipient'}, inplace=True)

result_anno = pd.DataFrame(columns=['id', 'sample', 'recipient_HGTC_n', 'recipient_HGTC_list', 'donor_HGTC_n', 'donor_HGTC_list', 'recipient', 'insert_locus', 'donor', 'delete_start', 'delete_end', 'reverse_flag'])
all_HGTC_idx = []
for idx in df.index:
    recipient_df, donor_df = search_row(idx, df, db, fr_size)
    all_HGTC_idx.extend(recipient_df.index)
    all_HGTC_idx.extend(donor_df.index)
    id = 'HGT_c{}'.format(idx+1)
    sample = df.loc[idx, 'sample']
    recipient_HGTC_n = recipient_df.shape[0]
    recipient_HGTC_list = ';'.join(recipient_df.index)
    if recipient_HGTC_n == 0:
        recipient_HGTC_list = 'NA'
    donor_HGTC_n = donor_df.shape[0]
    donor_HGTC_list = ';'.join(donor_df.index)
    if donor_HGTC_n == 0:
        donor_HGTC_list = 'NA'
    recipient = df.loc[idx, 'recipient']
    insert_locus = df.loc[idx, 'insert_locus']
    donor = df.loc[idx, 'donor']
    delete_start = df.loc[idx, 'delete_start']
    delete_end = df.loc[idx, 'delete_end']
    reverse_flag = df.loc[idx, 'reverse_flag']
    result_anno.loc[len(result_anno), ] = [id, sample, recipient_HGTC_n, recipient_HGTC_list, donor_HGTC_n, donor_HGTC_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]
    #result_anno.iloc[len(result_anno), ] = [id, sample, recipient_HGTC_n, recipient_HGTC_list, donor_HGTC_n, donor_HGTC_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]

db.loc[all_HGTC_idx, ].to_csv(os.path.join(outdir, 'output.HGTC_annotation.HGTC.tsv'), index=True, sep='\t')
result_anno.to_csv(os.path.join(outdir, 'output.HGTC_annotation.annotated.tsv'), index=False, sep='\t')
