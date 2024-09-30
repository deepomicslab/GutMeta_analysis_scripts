import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import networkx as nx
import pandas as pd
import pickle
# from pyfaidx import Fasta
# from sklearn.cluster import DBSCAN
from collections import Counter, defaultdict
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
# from Bio import SeqIO
from collections import Counter

COG_annotation = ["A: RNA processing and modification", "B: Chromatin structure and dynamics", "C: Energy production and conversion", "D: Cell cycle control, cell division, chromosome partitioning", "E: Amino acid transport and metabolism", "F: Nucleotide transport and metabolism", "G: Carbohydrate transport and metabolism", "H: Coenzyme transport and metabolism", "I: Lipid transport and metabolism", "J: Translation, ribosomal structure and biogenesis", "K: Transcription", "L: Replication, recombination and repair", "M: Cell wall/membrane/envelope biogenesis", "N: Cell motility", "O: Posttranslational modification, protein turnover, chaperones", "P: Inorganic ion transport and metabolism", "Q: Secondary metabolites biosynthesis, transport and catabolism", "R: General function prediction only", "S: Function unknown", "T: Signal transduction mechanisms", "U: Intracellular trafficking, secretion, and vesicular transport", "V: Defense mechanisms", "W: Extracellular structures", "Y: Nuclear structure", "Z: Cytoskeleton"]

def get_COG_dict():
    COG_dict = {}
    for anno in COG_annotation:
        arr = anno.split(":")
        name = arr[0]
        COG_dict[name] = anno  #arr[1].strip()

    COG_profile_dict = {}
    raw_dict = {"Metabolism":"QPIHFEGC", "Cellular processes and signaling":"XOUZNMTVDWY", "Information storage and Processing":"ABLKJ", "Function unknown":"RS"}
    for key in raw_dict:
        for i in range(len(raw_dict[key])):
            COG_profile_dict[raw_dict[key][i]] = key

    return COG_dict, COG_profile_dict

def enrichment_analysis(my_list, background_list, locus_type, data):
    COG_dict, COG_profile_dict = get_COG_dict()
    my_dict = Counter(my_list)
    background_dict = Counter(background_list)
    # data = []
    # for category in my_dict:
    for category in COG_dict:
        if category == "R" or category == "Y":
            continue
        if category in my_dict:
            a = my_dict[category]
        else:
            a = 0
        b = len(my_list) - a
        if category in background_dict:
            c = background_dict[category]
        else:
            c = 0
        d = len(background_list) - c

        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
    #    print (category, p_value, oddsratio, a, b, c, d) 
        data.append([category, COG_dict[category], p_value, oddsratio, a, COG_profile_dict[category], locus_type])
    return data


def cog_enrich(opath, group1_cogs, group2_cogs):
        ## COG enrichment analysis
    data = enrichment_analysis(group1_cogs, group2_cogs, "cog", [])
    df = pd.DataFrame(data, columns = ["category", "category_detail", "p_value", "fold", "gene_num", "profile", "locus_type"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected
    df.to_csv(opath, sep='\t', index=False)
    # print ("enriched COG num", len(data))