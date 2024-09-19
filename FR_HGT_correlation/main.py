import os
import sys  
import pandas as pd
import hgt_gcn
import util
import copy
import getopt


'''
    This is nfr adjuested by HGT comparison.
    options:
    --abdf   <str> input file path of related abundance
    --gcn_d    <str> input GCN distance
    --gcn    <str> input GCN
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --odir    <str> output directory
    --sp_gcn <str> input file of species genome annotation
    --ann    <str> input file of group info
    --groupid  <str> column name used for grouping, default: phenotype
    --method <str> method for correlation [pearson/spearman]

'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

odir = '.'
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'gcn_d=', 'gcn=', 'top_n=', 'db_dir=', 'hgt=', 'odir=', 'sp_g=', 'ann=', 'groupid=', 'method='])
for op, arg in ops:
    if op == '--abdf':
        abd_file = arg
    if op == '--gcn_d':
        d_file = arg
    if op == '--gcn':
        gcn_file = arg
    if op == '--odir':
        odir = arg
    # if op == '--top_n':
    #     top_n = float(arg)
    if op == '--db_dir':
        db_dir = arg
    if op == '--sp_g':
        sp_file = arg
    if op == '--hgt':
        hgt_file = arg
    if op == '--ann':
        metadata_file = arg
    if op == '--groupid':
        groupid = arg
    if op == '--method':
        method = arg


if not os.path.exists(odir):
    os.makedirs(odir)
top_n = 0
hgt_df = pd.read_csv(hgt_file, index_col=None, header=0)
group = pd.read_csv(metadata_file, index_col=0, header=0, sep='\t')
abd_df = pd.read_csv(abd_file, index_col=0, header=0, sep='\t')
sp_df = pd.read_csv(sp_file, index_col=0, header=0, sep='\t')
gcn_df = pd.read_csv(gcn_file, index_col=0, header=0, sep='\t')
sp_d = pd.read_csv(d_file, index_col=0, header=0, sep='\t')

group = group[[groupid]]
if not util.check_valid(group, abd_df):
    exit(2)

pheno_set = list(set(group[groupid]))
pheno_samples = {}
for g in pheno_set:
    pheno_samples[g] = list(group[group[groupid] == g].index)

abd_df = hgt_gcn.multi_sample_normalize(abd_df)
genome_ko = hgt_gcn.ko_df(hgt_df, db_dir)
sp_ko_df = hgt_gcn.hgt2sp_ko(sp_df, genome_ko)
hgt_nets = hgt_gcn.hgt2sp_hgt(sp_df, genome_ko)


nfr_result_df = pd.DataFrame(columns=['sample', 'nFR', 'adj_nFR', 'group'])
sum_fr_dict = {}
sum_adj_fr_dict = {}
sum_hgt_net_dict ={}

for g, slist in pheno_samples.items():
    # multi sample test
    sum_fr_net = pd.DataFrame()
    sum_adj_fr_net = pd.DataFrame()
    sum_hgt_net_df = pd.DataFrame()
    for sname in slist:
        if sname in hgt_nets.keys():
            sum_hgt_net_df = hgt_gcn.net_sum(sum_hgt_net_df, hgt_nets[sname])
        nfr_result_df.loc[sname, 'sample'] = sname
        part_abd_df = abd_df[sname]
        part_abd_df = part_abd_df[part_abd_df > 0]
        tmp_abd = list(part_abd_df.index)
        part_df = sp_ko_df[sp_ko_df['sample'] == sname][['sp1', 'sp2', 'ko', 'num']]
        common_sp = list(set(tmp_abd).intersection(set(gcn_df.index)))
        tmp_d = sp_d.loc[common_sp, common_sp]
        # original fr
        nfr_value, fr_df, profile = hgt_gcn.nfr(tmp_d, abd_df, sname)
        nfr_result_df.loc[sname, 'nFR'] = nfr_value
        nfr_result_df.loc[sname, 'group'] = g
        # align and add to sum nfr net
        sum_fr_net = hgt_gcn.net_sum(sum_fr_net, fr_df)
        if len(part_df)>0:
            new_gcn_df, effect_list = hgt_gcn.hgt_adjust_gcn(gcn_df, part_df)
            effect_list = list(set(effect_list).intersection(set(common_sp)))
            tmp_gcn = gcn_df.T[common_sp]
            if len(effect_list) < 20:
                new_d = hgt_gcn.adjust_d(tmp_d, tmp_gcn, effect_list)
            else:
                new_d = hgt_gcn.make_d(new_gcn_df.loc[common_sp,])
        
            nfr_value, fr_df, profile = hgt_gcn.nfr(new_d, abd_df, sname)
            nfr_result_df.loc[sname, 'adj_nFR'] = nfr_value
            sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
        else:
            nfr_result_df.loc[sname, 'adj_nFR'] = nfr_value
            sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
    sum_fr_dict[g] = copy.deepcopy(sum_fr_net)
    sum_adj_fr_dict[g] = copy.deepcopy(sum_adj_fr_net)
    sum_hgt_net_dict[g] = copy.deepcopy(sum_hgt_net_df)
    common_sp = list(set(sum_fr_net.index).intersection(set(sum_hgt_net_df.index)))
    if len(common_sp) == 0:
        print("No HGT found for current speceis")
    sum_hgt_net_df = sum_hgt_net_df.loc[common_sp, common_sp]
    mask = copy.deepcopy(sum_hgt_net_df)
    mask[mask > 0] = 1
    sum_fr_net = sum_fr_net.loc[common_sp, common_sp].multiply(mask)
    sum_adj_fr_net = sum_adj_fr_net.loc[common_sp, common_sp].multiply(mask)
    
    output1 = os.path.join(odir, 'output.fr_hgt_corr.sum_nFR.{}.tsv'.format(g))
    output2 = os.path.join(odir, 'output.fr_hgt_corr.sum_adj_nFR.{}.tsv'.format(g))
    output3 = os.path.join(odir, 'output.fr_hgt_corr.sum_hgt.{}.tsv'.format(g))
    output_fr_net = hgt_gcn.output_fr_net(sum_fr_net, top_n)[0]
    output_adj_fr_net = hgt_gcn.output_fr_net(sum_adj_fr_net, top_n)[0]
    output_hgt_net = hgt_gcn.output_fr_net(sum_hgt_net_df, top_n)[0]
    output_fr_net.columns = ['species1', 'species2', 'weight']
    output_adj_fr_net.columns = ['species1', 'species2', 'weight']
    output_hgt_net.columns = ['species1', 'species2', 'weight']
    output_fr_net.to_csv(output1, sep='\t', index=False)
    output_adj_fr_net.to_csv(output2, sep='\t', index=False)
    output_hgt_net.to_csv(output3, sep='\t', index=False)


result_df = pd.DataFrame(columns=['group', 'nFR', 'adj_nFR'])
for g in pheno_set:
    result_df.loc[g, 'group'] = g
    result_df.loc[g, 'nFR-HGT'] = hgt_gcn.net_correlation(sum_fr_dict[g], sum_hgt_net_dict[g], method)
    result_df.loc[g, 'adj_nFR-HGT'] = hgt_gcn.net_correlation(sum_adj_fr_dict[g], sum_hgt_net_dict[g], method)

output4 = os.path.join(odir, 'output.fr_hgt_corr.correlation.tsv')
result_df.to_csv(output4, sep='\t', index=False)