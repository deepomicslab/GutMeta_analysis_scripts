import gzip
import os
import pandas as pd
from collections import Counter
import copy
import numpy as np
from scipy.spatial.distance import pdist, squareform

# hgt_kegg_df: pandas.DataFrame sp1, sp2, ko, num (sorted sp1, sp2)
# gcn_df: pandas.DataFrame column is kos, index are species
def hgt_adjust_gcn(gcn_df, hgt_kegg_df):
    for idx in hgt_kegg_df.index:
        sp1, sp2, ko, num = hgt_kegg_df.loc[idx]
        if sp1 not in gcn_df.index or sp2 not in gcn_df.index or ko not in gcn_df.columns:
            continue
        gcn_df.loc[sp1, ko] += num
        gcn_df.loc[sp2, ko] += num
    return gcn_df

# extract ko
def gz2df(ifile):
    with gzip.open(ifile, 'rb') as f:
        df = pd.read_csv(f, header=None, index_col=None, sep='\t')
    return df

def overlap(range1, range2):
    if range1[0] > range2[1] or range1[1] < range2[0]:
        return False
    else:
        return True
    
def search_event(scaffold, db, range):
    scaffold_df = db[db[0] == scaffold]
    valid_idx = []
    for idx in scaffold_df.index:
        if overlap(range, [scaffold_df.loc[idx, 3], scaffold_df.loc[idx, 4]]):
            valid_idx.append(idx)
    valid_df = scaffold_df.loc[valid_idx, ]
    return valid_df

def search_row(idx, df, db_dir):
    tmp = '{}.gff.gz'
    row = df.loc[idx, ]
    # for donor
    donor = row['donor']
    range = [row['delete_start'] , row['delete_end']]
    chrom = donor.split('_')[0]
    ifile = os.path.join(db_dir, tmp.format(chrom))
    if not os.path.exists(ifile):
        return None
    db = gz2df(ifile)
    donor_df = search_event(donor, db, range)
    return donor_df

def extract_line(content):
    ko_list = []
    # check kegg
    if 'KEGG=' in content:
        kegg = content.split('KEGG=')[1].strip().split(';')[0]
        kos = kegg.split(',')
        #print('kos', kos)
        for ko in kos:
            if ko.startswith('ko:'):
                ko = ko.split('ko:')[-1]
                ko_list.append(ko)
    return ko_list

def extract(df):
    ko_list = []
    for idx in df.index:
        kos = extract_line(df.loc[idx, 8])
        ko_list += kos
    return ko_list

def ko_df(df, db_idir):
    df.rename(columns={'receptor':'recipient'}, inplace=True)
    result_anno = pd.DataFrame(columns=['id', 'sample', 
                                    'donor_KEGG_n', 'donor_KEGG_list', 
                                    'donor', 'recipient'])
    for idx in df.index:
        donor_df = search_row(idx, df, db_idir)
        if donor_df is None:
            continue
        id = 'HGT_c{}'.format(idx+1)
        sample = df.loc[idx, 'sample']
        ko_list = extract(donor_df)
        donor_KEGG_n = len(set(ko_list))
        donor_KEGG_list = ';'.join(set(ko_list))
        if donor_KEGG_n == 0:
            donor_KEGG_list = 'NA'
        donor = df.loc[idx, 'donor']
        recipient = df.loc[idx, 'recipient']
        result_anno.loc[len(result_anno), ] = [id, sample, donor_KEGG_n, donor_KEGG_list, donor, recipient]
    return result_anno

# sp_df : genome_id, species
def hgt2sp_ko(sp_df, result_anno):
    sp_ko_df = pd.DataFrame(columns=['sample', 'sp1', 'sp2', 'ko', 'num'])
    result_dict = {}
    for idx in result_anno.index:
        donor = result_anno.loc[idx, 'donor'].split('_')[0].split('.')[0]
        recipient = result_anno.loc[idx, 'recipient'].split('_')[0].split('.')[0]
        sample = result_anno.loc[idx, 'sample']
        if sample not in result_dict.keys():
            result_dict[sample] = {}
        if donor not in sp_df.index or recipient not in sp_df.index:
            continue
        sp1 = sp_df.loc[donor, 'species']
        sp2 = sp_df.loc[recipient, 'species']
        sp_pair = sorted([sp1, sp2])
        sp1 = sp_pair[0]
        sp2 = sp_pair[1]
        if sp1 == sp2:
            continue
        if sp1 not in result_dict[sample].keys():
            result_dict[sample][sp1] = {}
        if sp2 not in result_dict[sample][sp1].keys():
            result_dict[sample][sp1][sp2] = []
        ko_list = result_anno.loc[idx, 'donor_KEGG_list'].split(';')
        result_dict[sample][sp1][sp2] += ko_list
    for sample in result_dict.keys():
        for sp1 in result_dict[sample].keys():
            for sp2 in result_dict[sample][sp1].keys():
                ko_list = result_dict[sample][sp1][sp2]
                ko_count = Counter(ko_list)
                for ko in ko_count.keys():
                    sp_ko_df.loc[len(sp_ko_df), ] = [sample, sp1, sp2, ko, ko_count[ko]]
  
    return sp_ko_df


# extract HGT species network
def hgt2sp_hgt(sp_df, result_anno):
    hgt_dict = {}
    result_dict = {}
    for idx in result_anno.index:
        donor = result_anno.loc[idx, 'donor'].split('_')[0].split('.')[0]
        recipient = result_anno.loc[idx, 'recipient'].split('_')[0].split('.')[0]
        sample = result_anno.loc[idx, 'sample']
        if sample not in hgt_dict.keys():
            hgt_dict[sample] = {}
        if donor not in sp_df.index or recipient not in sp_df.index:
            continue
        sp1 = sp_df.loc[donor, 'species']
        sp2 = sp_df.loc[recipient, 'species']
        if sp1 == sp2:
            continue
        sp_pair = sorted([sp1, sp2])
        sp1 = sp_pair[0]
        sp2 = sp_pair[1]
        if sp1 not in hgt_dict[sample].keys():
            hgt_dict[sample][sp1] = {}
        if sp2 not in hgt_dict[sample][sp1].keys():
            hgt_dict[sample][sp1][sp2] = 0
        hgt_dict[sample][sp1][sp2] += 1
       
    for sample in hgt_dict.keys():
        new_df = pd.DataFrame()
        for sp1 in hgt_dict[sample].keys():
            for sp2 in hgt_dict[sample][sp1].keys():
                new_df.loc[sp1, sp2] = hgt_dict[sample][sp1][sp2]
                new_df.loc[sp2, sp1] = hgt_dict[sample][sp1][sp2]
        new_df.fillna(0, inplace=True)
        result_dict[sample] = copy.deepcopy(new_df)
    return result_dict
        

def nfr(d_df, profile, sname):
    sp_list = list(set(list(profile.index)).intersection(set(list(d_df.index))))
    n = len(sp_list)
    corr = np.ones(shape=(n, n)) - d_df.loc[sp_list, sp_list].values
    np.fill_diagonal(corr, 0)
    a = np.array(profile.loc[sp_list, sname])
    inter_matrix = np.dot(a.reshape(len(a), 1),a.reshape(1, len(a)))
    np.fill_diagonal(inter_matrix, 0)
    td = np.sum(inter_matrix)/2
    fr = np.sum(np.multiply(inter_matrix, corr))/2
    fr_df = pd.DataFrame(np.multiply(inter_matrix, corr), index=sp_list, columns=sp_list)
    profile = profile.loc[sp_list, sname]
    if td == 0:
        return 0
    return fr/td, fr_df, profile

def Jaccard_distance(x, y):
    union = sum(np.maximum(x, y))
    intersection = sum(np.minimum(x, y))
    # print(intersection, union)
    if union != 0:
        return 1 - intersection/union
    else:
        # print("union size is 0")
        return 1

def make_d(ref_GCN, sep='\t'):
    ref_GCN = ref_GCN.T
    sp_list = list(ref_GCN.columns)
    distance_compressed = pdist(ref_GCN.values.T, Jaccard_distance)
    distance_matrix = squareform(distance_compressed)
    distance_df = pd.DataFrame(distance_matrix, columns = sp_list, index = sp_list)
    return distance_df            
        
