import sys
import subprocess
import getopt
import os
import pandas as pd
import gzip
import COG_enrichment as ce
import kegg_enrichment as ke
import pickle

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
    recipient_excluded_df = get_backgroud(recipient, db, recipient_range)
    # for donor
    donor = row['donor']
    range = [max(0, row['delete_start'] - fr_size), row['delete_end'] + fr_size]
    chrom = donor.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    db = gz2df(ifile)
    donor_df = search_event(donor, db, range)
    donor_excluded_df = get_backgroud(donor, db, range)
    return recipient_df, donor_df, recipient_excluded_df, donor_excluded_df 

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
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --fr_size    <int> flanking region size
    --ko_pathway_dict <str> input file of ko_pathway_dict
    --outdir <str> output dir
'''

ops, args = getopt.getopt(sys.argv[1:], '', ['db_dir=', 'hgt=', 'fr_size=', 'ko_pathway_dict=', 'outdir='])
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

with open(pfile, 'rb') as f: 
    ko_pathway_dict = pickle.load(f)
df = pd.read_csv(infile, header=0, index_col=None)
df.rename(columns={'receptor':'recipient'}, inplace=True)

result_anno = pd.DataFrame(columns=['id', 'sample', 'recipient_KEGG_n', 'recipient_KEGG_list', 'recipient_COG_n', 'recipient_COG_list',
                                    'donor_KEGG_n', 'donor_KEGG_list', 'donor_COG_n', 'donor_COG_list',
                                    'recipient', 'insert_locus', 'donor', 'delete_start', 'delete_end', 'reverse_flag'])

bk_ko_set = set()
bk_cog_set = set()
exist_ko_set = set()
exist_cog_set = set()


for idx in df.index:
    recipient_df, donor_df, recipient_excluded_df, donor_excluded_df = search_row(idx, df, db_idir, fr_size)
    id = 'HGT_c{}'.format(idx+1)
    sample = df.loc[idx, 'sample']
    ko_list, cog_list, cog_cate = extract(recipient_df)
    exist_ko_set.update(set(ko_list))
    exist_cog_set.update(set(cog_cate))
    recipient_KEGG_n = len(set(ko_list))
    recipient_COG_n = len(set(cog_list))
    recipient_KEGG_list = ';'.join(set(ko_list))
    recipient_COG_list = ';'.join(set(cog_list))
    if recipient_KEGG_n == 0:
        recipient_KEGG_list = 'NA'
    if recipient_COG_n == 0:
        recipient_COG_list = 'NA'

    ko_list, cog_list, cog_cate = extract(donor_df)
    exist_ko_set.update(set(ko_list))
    exist_cog_set.update(set(cog_cate))
    donor_KEGG_n = len(set(ko_list))
    donor_COG_n = len(set(cog_list))
    donor_KEGG_list = ';'.join(set(ko_list))
    donor_COG_list = ';'.join(set(cog_list))
    if donor_KEGG_n == 0:
        donor_KEGG_list = 'NA'
    if donor_COG_n == 0:
        donor_COG_list = 'NA'
    recipient = df.loc[idx, 'recipient']
    insert_locus = df.loc[idx, 'insert_locus']
    donor = df.loc[idx, 'donor']
    delete_start = df.loc[idx, 'delete_start']
    delete_end = df.loc[idx, 'delete_end']
    reverse_flag = df.loc[idx, 'reverse_flag']
    result_anno.loc[len(result_anno), ] = [id, sample, recipient_KEGG_n, recipient_KEGG_list, recipient_COG_n, recipient_COG_list, donor_KEGG_n, donor_KEGG_list, donor_COG_n, donor_COG_list, recipient, insert_locus, donor, delete_start, delete_end, reverse_flag]

    # add bg ko
    ko_list, cog_list, cog_cate = extract(recipient_excluded_df)
    bk_ko_set.update(set(ko_list))
    bk_cog_set.update(set(cog_cate))

    ko_list, cog_list, cog_cate = extract(donor_excluded_df)
    bk_ko_set.update(set(ko_list))
    bk_cog_set.update(set(cog_cate))
    
result_anno.to_csv(os.path.join(outdir, 'output.functional_annotation.annotated.tsv'), index=False, sep='\t')
background_counts = ke.get_pathways(list(bk_ko_set), ko_pathway_dict)
input_counts = ke.get_pathways(list(exist_ko_set), ko_pathway_dict)
opath = os.path.join(outdir, 'output.functional_annotation.enrich.KEGG.tsv')
ke.enrichment_analysis(list(exist_ko_set), list(bk_ko_set), input_counts, background_counts, opath)
opath = os.path.join(outdir, 'output.functional_annotation.enrich.COG.tsv')
ce.cog_enrich(opath, list(exist_cog_set), list(bk_cog_set) )