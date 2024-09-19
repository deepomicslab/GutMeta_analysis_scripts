# compute old FR average and new FR average network and nFR value 
# param nedge
import os
import sys  
import pandas as pd
import hgt_gcn
import getopt

'''
    This is to compute new nfr adjusted by HGT.
    options:
    --abdf   <str> input file path of related abundance
    --gcn_d    <str> input GCN distance
    --gcn    <str> input GCN
    --top_n    <str> output top n edges
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --odir    <str> output directory
    --sp_g   <str> input file of species genome annotation
'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

odir = '.'
top_n = 0
ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'gcn_d=', 'gcn=', 'top_n=', 'db_dir=', 'hgt=', 'odir=', 'sp_g='])
for op, arg in ops:
    if op == '--abdf':
        abd_file = arg
    if op == '--gcn_d':
        d_file = arg
    if op == '--gcn':
        gcn_file = arg
    if op == '--odir':
        odir = arg
    if op == '--top_n':
        top_n = float(arg)
    if op == '--db_dir':
        db_dir = arg
    if op == '--sp_g':
        sp_file = arg
    if op == '--hgt':
        hgt_file = arg


if not os.path.exists(odir):
    os.makedirs(odir)

hgt_df = pd.read_csv(hgt_file, index_col=None, header=0)
abd_df = pd.read_csv(abd_file, index_col=0, header=0, sep='\t')
sp_df = pd.read_csv(sp_file, index_col=0, header=0, sep='\t')
gcn_df = pd.read_csv(gcn_file, index_col=0, header=0, sep='\t')
sp_d = pd.read_csv(d_file, index_col=0, header=0, sep='\t')

abd_df = hgt_gcn.multi_sample_normalize(abd_df)
genome_ko = hgt_gcn.ko_df(hgt_df, db_dir)
sp_ko_df = hgt_gcn.hgt2sp_ko(sp_df, genome_ko)

# multi sample adjust
nfr_result_df = pd.DataFrame(columns=['sample', 'nFR', 'adj_nFR'])
sum_fr_net = pd.DataFrame()
sum_adj_fr_net = pd.DataFrame()
for sname in list(abd_df.columns):
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
    # align and add to sum nfr net
    sum_fr_net = hgt_gcn.net_sum(sum_fr_net, fr_df)
    if len(part_df)>0:
        new_gcn_df, effect_list = hgt_gcn.hgt_adjust_gcn(gcn_df, part_df)
        effect_list = list(set(effect_list).intersection(set(common_sp)))
        tmp_gcn = gcn_df.T[common_sp]
        if len(effect_list) < 1000:
            new_d = hgt_gcn.adjust_d(tmp_d, tmp_gcn, effect_list)
        else:
            new_d = hgt_gcn.make_d(new_gcn_df.loc[common_sp,])
    
        nfr_value, fr_df, profile = hgt_gcn.nfr(new_d, abd_df, sname)
        nfr_result_df.loc[sname, 'adj_nFR'] = nfr_value
        sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
    else:
        nfr_result_df.loc[sname, 'adj_nFR'] = nfr_value
        sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
avg_fr_net = sum_fr_net / len(abd_df.columns)
avg_adj_fr_net = sum_adj_fr_net / len(abd_df.columns)

# top n 
avg_fr_output = hgt_gcn.output_fr_net(avg_fr_net, top_n)[0]
avg_adj_fr_output = hgt_gcn.output_fr_net(avg_adj_fr_net, top_n)[0]

output_path1 = os.path.join(odir, 'output.hgt_adjusted_fr.avg_nFR.tsv')
output_path2 = os.path.join(odir, 'output.hgt_adjusted_fr.avg_adj_nFR.tsv')
output_path3 = os.path.join(odir, 'output.hgt_adjusted_fr.nFR_adjusted.tsv')

nfr_result_df.to_csv(output_path3, sep='\t', index=False)
avg_fr_output.columns = ['species1', 'species2', 'weight']
avg_adj_fr_output.columns = ['species1', 'species2', 'weight']
avg_fr_output.to_csv(output_path1, sep='\t', index=False)
avg_adj_fr_output.to_csv(output_path2, sep='\t', index=False)