# compute old FR average and new FR average network and nFR value 
# param nedge
import os
import sys  
import pandas as pd
import hgt_gcn
import util
import copy
import getopt
from util import *


#python main.py --abdf ../../hgt_abd_metadata/ERP010700.merged.tsv --ann ../../hgt_abd_metadata/ERP010700.metadata.v2.tsv --gcn_d ../../sp_d.tsv --groupid phenotype --method wilcox.test --odir .
'''
    This is nfr adjuested by HGT comparison.
    options:
    --abdf   <str> input file path of related abundance
    --gcn_d    <str> input GCN distance
    --gcn    <str> input GCN
    --top_n    <str> output top n edges
    --db_dir   <str> input dir of MGE database
    --hgt    <str> input file of HGT output
    --odir    <str> output directory
    --sp_gcn <str> input file of species genome annotation
    --ann    <str> input file of group info
    --groupid  <str> column name used for grouping, default: phenotype
    --method <str> method for testing [wilcox.test/t.test/kruskal.test/aov]

'''

python_file = os.path.abspath(__file__) 
python_dir = os.path.dirname(python_file)

odir = '.'
top_n = 0
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
    if op == '--top_n':
        top_n = float(arg)
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

hgt_df = pd.read_csv(hgt_file, index_col=None, header=0)
metadata_df = pd.read_csv(metadata_file, index_col=None, header=0, sep='\t')
abd_df = pd.read_csv(abd_file, index_col=0, header=0, sep='\t')
sp_df = pd.read_csv(sp_file, index_col=0, header=0, sep='\t')
gcn_df = pd.read_csv(gcn_file, index_col=0, header=0, sep='\t')
sp_d = pd.read_csv(d_file, index_col=0, header=0, sep='\t')

group = util.metadata2gf(metadata_file,groupid)
if not util.check_valid(group, abd_df):
    exit(2)

pheno_set = list(set(group[groupid]))
g1 = pheno_set[0]
g2 = pheno_set[1]
pheno_samples = {}
pheno_samples[g1] = list(group[group[groupid] == g1].index)
pheno_samples[g2] = list(group[group[groupid] == g2].index)

abd_df = hgt_gcn.multi_sample_normalize(abd_df)
genome_ko = hgt_gcn.ko_df(hgt_df, db_dir)
sp_ko_df = hgt_gcn.hgt2sp_ko(sp_df, genome_ko)

nfr_result_df = pd.DataFrame(columns=['sample', 'group', 'nFR', 'aFR'])
avg_fr_dict = {}
avg_adj_fr_dict = {}
for g, slist in pheno_samples.items():
    # multi sample test
    sum_fr_net = pd.DataFrame()
    sum_adj_fr_net = pd.DataFrame()
    for sname in slist:
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
            nfr_result_df.loc[sname, 'aFR'] = nfr_value
            sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
        else:
            nfr_result_df.loc[sname, 'aFR'] = nfr_value
            sum_adj_fr_net = hgt_gcn.net_sum(sum_adj_fr_net, fr_df)
    avg_fr_net = sum_fr_net / len(slist)
    avg_adj_fr_net = sum_adj_fr_net / len(slist)

    avg_fr_dict[g] = copy.deepcopy(avg_fr_net)
    avg_adj_fr_dict[g] = copy.deepcopy(avg_adj_fr_net)
    avg_fr_output = hgt_gcn.output_fr_net(avg_fr_net, top_n)[0]
    avg_adj_fr_output = hgt_gcn.output_fr_net(avg_adj_fr_net, top_n)[0]
    if g == g1:
        group_str = 'group1'
    else:
        group_str = 'group2'
    output1 = os.path.join(odir, 'output.aFR_comparison.nFR_average_network.{}.tsv'.format(group_str))
    output2 = os.path.join(odir, 'output.aFR_comparison.aFR_average_network.{}.tsv'.format(group_str))

    output12 = os.path.join(odir, 'plot_FR_average_network.{}.json'.format(group_str))
    output22 = os.path.join(odir, 'plot_aFR_average_network.{}.json'.format(group_str))

    plot_output_aFR_nFR_average_network_tsv1 = read_tsv(output1)
    plot_output_aFR_nFR_average_network1 = {'plotType': 'FR', 'data':  plot_output_aFR_nFR_average_network_tsv1[1:], 'columns': plot_output_aFR_nFR_average_network_tsv1[0]}
    save_dict_as_json(plot_output_aFR_nFR_average_network1, output12)

    plot_output_aFR_nFR_average_network_tsv2 = read_tsv(output2)
    plot_output_aFR_nFR_average_network2 = {'plotType': 'FR', 'data':  plot_output_aFR_nFR_average_network_tsv2[1:], 'columns': plot_output_aFR_nFR_average_network_tsv2[0]}
    save_dict_as_json(plot_output_aFR_nFR_average_network2, output22)

    avg_fr_output.columns = ['species1', 'species2', 'weight']
    avg_adj_fr_output.columns = ['species1', 'species2', 'weight']
    avg_fr_output.to_csv(output1, sep='\t', index=False)
    avg_adj_fr_output.to_csv(output2, sep='\t', index=False)


cols = [
    'group1',
    'group2',
    'g1_mean',
    'g2_mean',
    'g1_variance',
    'g2_variance',
    'g1_n',
    'g2_n',
    'p_value',
    'g1/g2', 
    'enriched']


p_df = pd.DataFrame(columns=cols)
p_df.loc['nFR', 'group1'] = g1
p_df.loc['nFR', 'group2'] = g2
g1_v = nfr_result_df[nfr_result_df['group'] == g1]['nFR'].values.astype(float)
g2_v = nfr_result_df[nfr_result_df['group'] == g2]['nFR'].values.astype(float)
p_df.loc['nFR', 'g1_mean'] = g1_v.mean()
p_df.loc['nFR', 'g2_mean'] = g2_v.mean()
if p_df.loc['nFR', 'g1_mean'] > p_df.loc['nFR', 'g2_mean']:
    p_df.loc['nFR', 'enriched'] = g1
else:
    p_df.loc['nFR', 'enriched'] = g2
p_df.loc['nFR', 'g1/g2'] = p_df.loc['nFR', 'g1_mean']/p_df.loc['nFR', 'g2_mean']
p_df.loc['nFR', 'g1_variance'] = g1_v.var()
p_df.loc['nFR', 'g2_variance'] = g2_v.var()
# p_df.loc[0, 'g1_occ'] = len(g1_v[g1_v > 0])/len(g1_v)
# p_df.loc[0, 'g2_occ'] = len(g2_v[g2_v > 0])/len(g2_v)
p_df.loc['nFR', 'g1_n'] = len(g1_v)
p_df.loc['nFR', 'g2_n'] = len(g2_v)
p_df.loc['nFR', 'p_value'] = hgt_gcn.test(method, g1_v, g2_v)


p_df.loc['aFR', 'group1'] = g1
p_df.loc['aFR', 'group2'] = g2
g1_v = nfr_result_df[nfr_result_df['group'] == g1]['aFR'].values.astype(float)
g2_v = nfr_result_df[nfr_result_df['group'] == g2]['aFR'].values.astype(float)
p_df.loc['aFR', 'g1_mean'] = g1_v.mean()
p_df.loc['aFR', 'g2_mean'] = g2_v.mean()
if p_df.loc['aFR', 'g1_mean'] > p_df.loc['aFR', 'g2_mean']:
    p_df.loc['aFR', 'enriched'] = g1
else:
    p_df.loc['aFR', 'enriched'] = g2
p_df.loc['aFR', 'g1/g2'] = p_df.loc['aFR', 'g1_mean']/p_df.loc['aFR', 'g2_mean']
p_df.loc['aFR', 'g1_variance'] = g1_v.var()
p_df.loc['aFR', 'g2_variance'] = g2_v.var()
#p_df.loc[1, 'g1_occ'] = len(g1_v[g1_v > 0])/len(g1_v)
#p_df.loc[1, 'g2_occ'] = len(g2_v[g2_v > 0])/len(g2_v)
p_df.loc['aFR', 'g1_n'] = len(g1_v)
p_df.loc['aFR', 'g2_n'] = len(g2_v)
p_df.loc['aFR', 'p_value'] = hgt_gcn.test(method, g1_v, g2_v)
outpath3 = os.path.join(odir, 'output.aFR_comparison.pvalue.tsv')
outpath4 = os.path.join(odir, 'output.aFR_comparison.results.tsv')
p_df.to_csv(outpath3, sep='\t', index=True)
nfr_result_df.to_csv(outpath4, sep='\t', index=False)


# plot main data
plot_FR_adjusted_comparison_data_tsv = read_tsv(outpath4)
plot_FR_adjusted_comparison_info_tsv = read_tsv(outpath3)

outpath32 = os.path.join(odir, 'plot_FR_adjusted_comparison_info.js')
outpath42 = os.path.join(odir, 'plot_FR_adjusted_comparison_data.json')


plot_FR_adjusted_comparison_data = {'data':  plot_FR_adjusted_comparison_data_tsv[1:], 'columns':  plot_FR_adjusted_comparison_data_tsv[0]}
plot_FR_adjusted_comparison_info = {'data':  plot_FR_adjusted_comparison_info_tsv[1:], 'columns':  plot_FR_adjusted_comparison_info_tsv[0]}

save_dict_as_json(plot_FR_adjusted_comparison_data, outpath42)
save_dict_as_json(plot_FR_adjusted_comparison_info, outpath32)

