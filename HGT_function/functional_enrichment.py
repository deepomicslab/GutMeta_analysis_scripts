#!/usr/bin/env python3

"""
Analyze the function of HGT-related genes and analyze the HGT transfer patterns
"""

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
import sys
import argparse

import kegg_enrichment

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

class Acc_Bkp(object):

    def __init__(self, list):

        self.from_ref = list[0]
        self.from_bkp = int(list[1])
        self.from_side = list[2]
        self.from_strand = list[3]
        self.from_split_reads = int(list[12])

        self.to_ref = list[4]
        self.to_bkp = int(list[5])
        self.to_side = list[6]
        self.to_strand = list[7]
        self.to_split_reads = int(list[13])

        self.if_reverse = list[8]
        
        self.cross_split_reads = int(list[14])
        self.pair_end = int(list[15])

        
        self.from_ref_genome = "_".join(self.from_ref.split("_")[:-1])
        self.to_ref_genome = "_".join(self.to_ref.split("_")[:-1])
        self.abundance = None
        self.split_abundance = None

        self.read = None
        self.score = float(list[11])

        # self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        # self.pair_tag = "&".join(sorted([self.from_ref_genome, self.to_ref_genome]))

class Event(object):

    def __init__(self, array):
        self.sample = array[1]
        self.ins_genome = array[2]
        self.ins_genome_pure = "_".join(self.ins_genome.split("_")[:-1])
        self.ins_pos = int(array[3])
        self.del_genome = array[4]
        self.del_start = int(array[5])
        self.del_end = int(array[6])
        self.reverse_flag = array[7]
        self.IS_flag = None
        self.Transposon_flag = None

        bin_size = 100
        self.tag = self.ins_genome + "&" + str(round(self.ins_pos/bin_size)) + "&" + self.del_genome + "&" + \
            str(round(self.del_start/bin_size)) + "&" + str(round(self.ins_pos/bin_size))

    def check_IS(self, min_gene_frac, annotation):
        transfer_interval = [self.del_start, self.del_end]
        if self.del_genome in annotation.gene_annotation:
            self.IS_flag = True
            self.Transposon_flag = True
            # print ("Event", self.del_genome, transfer_interval)
            gene_intervals = annotation.gene_annotation[self.del_genome]["intervals"]
            for gene_interval in gene_intervals:
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], [transfer_interval], min_gene_frac)
                if not locate_transfer_flag:
                    continue
                gene_anno_dict = annotation.gene_annotation[self.del_genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "product" not in gene_anno_dict:
                    gene_anno_dict["product"] = "NA"
                if not re.search("IS[0-9]", gene_anno_dict["product"]):
                    self.IS_flag = False
                product_classification = annotation.classify_product(gene_anno_dict["product"])
                if product_classification != "Transposon":
                    self.Transposon_flag = False
                # print (gene_interval, gene_anno_dict["product"], self.IS_flag, self.Transposon_flag, sep = "\t")
        else:
            self.IS_flag = False
        # print ("Event", self.del_genome, transfer_interval, self.IS_flag)
        # print ("<<<<<<<<<<<<<<<<<<<<\n")

def get_gene_lengths(identified_hgt):
    HGT_event_list = []
    contig_lengths = []
    for line in open(identified_hgt):
        array = line.strip().split(",")
        if array[1] == "sample":
            continue
        gene_len = int(array[6]) - int(array[5])
        event = Event(array)
        HGT_event_list.append(event)
        # if gene_len == 1543502:
        #     print (array)
        contig_lengths.append(gene_len)

    print (np.mean(contig_lengths), np.median(contig_lengths), len(contig_lengths), min(contig_lengths), max(contig_lengths), sum(contig_lengths))
    print (sorted(contig_lengths)[-10:])

    data = []
    for a in contig_lengths:
        data.append([a, "mm"])
    df = pd.DataFrame(data, columns = ["Length", "mm"])
    df.to_csv("/mnt/d/R_script_files//gene_len.csv", sep=',')   

    print ("Finish reading events, num is", len(HGT_event_list))
    return HGT_event_list

class Transfer_times():

    def __init__(self):
        self.HGT_event_dict = {}
        self.max_diff = 50

    def read_events(self, identified_hgt):
        i = 0
        for line in open(identified_hgt):
            array = line.strip().split(",")
            array = ["x"] + array
            if array[1] == "sample":
                continue
            event = Event(array)
            sample = array[1]
            # if phenotype_dict[sample][0] == "Time-series":
            #     continue
            # print (phenotype_dict[sample][0])
            i += 1

            if sample not in self.HGT_event_dict:
                self.HGT_event_dict[sample] = []
            self.HGT_event_dict[sample].append(event) 
        print ("sample num ", len(self.HGT_event_dict), "kept HGTs num is", i)

    def read_events_bk(self, identified_hgt):
        i = 0
        for line in open(identified_hgt):
            array = line.strip().split(",")
            if array[1] == "sample":
                continue
            event = Event(array)
            sample = array[1]
            # if phenotype_dict[sample][0] == "Time-series":
            #     continue
            # print (phenotype_dict[sample][0])
            i += 1

            if sample not in self.HGT_event_dict:
                self.HGT_event_dict[sample] = []
            self.HGT_event_dict[sample].append(event) 
        print ("sample num ", len(self.HGT_event_dict), "kept HGTs num is", i)

class Annotation():

    def __init__(self, gff):
        self.gff = gff
        self.near = 100
        self.min_gene_frac = 0.5
        self.gene_annotation = {}
        self.gene_classification = {}
        self.pattern_dict = {}
        self.init_pattern_dict()

    def read_gff(self):

        f = open(self.gff)
        for line in f:
            array = line.split("\t")
            if len(array) < 9:
                continue
            genome = array[0]
            g_type = array[2]
            detail = array[8].strip()
            start = int(array[3])
            end = int(array[4])
            if genome not in self.gene_annotation:
                self.gene_annotation[genome] = {}
                self.gene_annotation[genome]["intervals"] = []
            self.gene_annotation[genome]["intervals"].append([start, end])
            self.gene_annotation[genome][str(start)+ "_"+str(end)] = self.understand_gene(detail)
        f.close() 

    def given_seg(self, genome, gene_interval):
        if genome not in self.gene_annotation:
            return ["NA"]
        intervals = self.gene_annotation[genome]["intervals"]
        genes_around = []
        for inter in intervals:
            gene_ID, product = '', ''
            gene_anno_dict = self.gene_annotation[genome][str(inter[0]) + "_" + str(inter[1])]
            if "ID" in gene_anno_dict:
                gene_ID = gene_anno_dict["ID"]
            if "product" in gene_anno_dict:
                product = gene_anno_dict["product"]

            CDS_length = int(inter[1]) - int(inter[0])
            if inter[0] >= gene_interval[0] and inter[0] <= gene_interval[1] and (gene_interval[1] - inter[0])/CDS_length > self.min_gene_frac :
                genes_around.append(product)
            elif inter[0] <= gene_interval[0] and inter[1] >= gene_interval[0] and (inter[1] - gene_interval[0])/CDS_length > self.min_gene_frac:
                genes_around.append(product)

        if len(genes_around) > 0:
            # print (genome, gene_interval, gene_interval[1] - gene_interval[0], CDS_length, genes_around)
            return genes_around
        else:
            return []

    def understand_gene(self, detail):
        array = detail.split(";")
        anno_dict = {}
        for arr in array:
            category = arr.split("=")
            anno_dict[category[0]] = category[1]
        return anno_dict

    def init_pattern_dict(self):
        Transposon = "transpos*; insertion; resolv*; Tra[A-Z]; Tra[0-9]; IS[0-9]; conjugate transposon"
        Plasmid = "resolv*; relax*; conjug*; trb; mob*; plasmid; type IV; toxin; chromosome partitioning; chromosome segregation"
        Phage ="capsid; phage; tail; head; tape measure; antitermination; antiterminatio"
        Other_HGT_machinery=" integrase; excision*; exonuclease; recomb; toxin; CRISPR; restrict*; resolv*; topoisomerase; reverse transcrip"
        
        Carbohydrate_active =  "Genes present in the CAZY database; glycosyltransferase; glycoside hydrolase; xylan; monooxygenase; rhamnos*; cellulose; sialidase; *ose; acetylglucosaminidase; cellobiose; galact*; fructose; aldose; starch; mannose; mannan*; glucan; lyase; glycosyltransferase; glycosidase; pectin; SusD; SusC; fructokinase; galacto*; arabino*"
        antibiotic_resistance =  "Genes present in the ARDB; multidrug; azole resistance; antibiotic resistance; TetR; tetracycline resistance; VanZ; betalactam*; beta-lactam; antimicrob*; lantibio*"
        # hypothetical_protein = "hypothetical protein"
        
        gene_classification = {}
        self.pattern_dict["Transposon"] = get(Transposon)
        self.pattern_dict["Plasmid"] = get(Plasmid)
        self.pattern_dict["Phage"] = get(Phage)
        self.pattern_dict["Other_HGT_mechanisms"] = get(Other_HGT_machinery)
        self.pattern_dict["CAZYmes"] = get(Carbohydrate_active)
        self.pattern_dict["Antibiotic resistance"] = get(antibiotic_resistance)
        # self.pattern_dict["hypothetical protein"] = get(hypothetical_protein)
    
    def classify_product(self, product):
        # Transposon pattern
        transposon_pattern = re.compile(r"transpos\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon")

        # Plasmid pattern
        plasmid_pattern = re.compile(r"relax\S*|conjug\S*|mob\S*|plasmid|type IV|chromosome partitioning|chromosome segregation")

        # Phage pattern
        phage_pattern = re.compile(r"capsid|phage|tail|head|tape measure|antiterminatio")

        # Other HGT mechanisms pattern
        hgt_pattern = re.compile(r"integrase|excision\S*|exonuclease|recomb|toxin|restrict\S*|resolv\S*|topoisomerase|reverse transcrip")

        # Carbohydrate active pattern
        carbohydrate_pattern = re.compile(r"glycosyltransferase|glycoside hydrolase|xylan|monooxygenase|rhamnos\S*|cellulose|sialidase|\S*ose($|\s|\-)|acetylglucosaminidase|cellobiose|galact\S*|fructose|aldose|starch|mannose|mannan\S*|glucan|lyase|glycosyltransferase|glycosidase|pectin|SusD|SusC|fructokinase|galacto\S*|arabino\S*")

        # Antibiotic resistance pattern
        antibiotic_pattern = re.compile(r"azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\S*|beta-lactam|antimicrob\S*|lantibio\S*")

        product_classification = 'unclassified'

        # Search for matches
        if plasmid_pattern.search(product):
            # print("Found a plasmid!")
            product_classification = "plasmid"
        if phage_pattern.search(product):
            # print("Found a phage!")
            product_classification = "phage"
        if transposon_pattern.search(product):
            # print("Found a transposon!")
            product_classification = "transposon"
        if hgt_pattern.search(product):
            # print("Found an HGT mechanism!")
            product_classification = "Other_HGT_mechanisms"
        if carbohydrate_pattern.search(product):
            # print("Found a carbohydrate active enzyme!")
            product_classification = "CAZYmes"
        if antibiotic_pattern.search(product):
            # print("Found an antibiotic resistance gene!")
            product_classification = "ARG"

        return product_classification
    
def get(x):
    x = x.split(";")
    # print (x)
    for i in range(len(x)):
        x[i] = x[i].strip()
        if x[i][0] == "“":
            x[i] = x[i][1:]
        if x[i][0] == "*":
            x[i] = x[i][1:]
        if x[i][-1] == "”":
            x[i] = x[i][:-1]
        if x[i][-1] == "*":
            x[i] = x[i][:-1]
    pat = ''
    for i in range(len(x)):
        pat += x[i] + "|"
    if pat[-1] == "|":
        pat = pat[:-1]
    return pat


def merge_intervals(intervals):
    # Sort the intervals by their start time
    intervals.sort(key=lambda x: x[0])
    
    # Initialize an empty list to store the merged intervals
    merged = []
    
    # Iterate through all the intervals
    for interval in intervals:
        # If the merged list is empty or the current interval doesn't overlap with the last merged interval
        if not merged or interval[0] > merged[-1][1]:
            merged.append(interval)
        else:
            # If the current interval overlaps with the last merged interval, merge them
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
    
    # Return the merged intervals
    return merged

class Extract_KO():

    def __init__(self, HGT_event_dict, surround=0):
        self.HGT_event_dict = HGT_event_dict
        self.min_gene_frac = min_gene_frac
        self.transfer_regions = {}
        self.insert_sites = {}
        # self.no_transfer_regions = {}
        # self.no_ins_regions= {}
        self.transfer_kos = []
        self.no_transfer_kos = []
        self.insert_kos = []
        self.no_insert_kos = []
        self.near = surround
    
    def classify_regions(self):
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:


                if event.del_genome not in self.transfer_regions:
                    self.transfer_regions[event.del_genome] = []
                if event.ins_genome not in self.insert_sites:
                    self.insert_sites[event.ins_genome] = []
                self.transfer_regions[event.del_genome].append([event.del_start, event.del_end])
                self.insert_sites[event.ins_genome].append(event.ins_pos)
        self.merge_transfer_region()
    
    def merge_transfer_region(self):
        for genome in self.transfer_regions:
            self.transfer_regions[genome] = merge_intervals(self.transfer_regions[genome])

    def classify_kos(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")
                new_KEGG_list = []
                for element in KEGG_list:
                    if element[:3] == "ko:":
                        element = element[3:]
                        new_KEGG_list.append(element)
                KEGG_list = new_KEGG_list
                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    self.transfer_kos += KEGG_list
                else:
                    self.no_transfer_kos += KEGG_list
                # print (KEGG_list, locate_transfer_flag)
            # break
        # print (len(self.transfer_kos), len(self.no_transfer_kos))
        return self.transfer_kos, self.no_transfer_kos

    def collect_all_ko(self):
        all_ko = []
        for genome in annotation.gene_annotation:
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")
                all_ko += KEGG_list
        print_data(all_ko, all_kos_file)

    def classify_kos_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")

                new_KEGG_list = []
                for element in KEGG_list:
                    if element[:3] == "ko:":
                        element = element[3:]
                        new_KEGG_list.append(element)
                KEGG_list = new_KEGG_list

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0]- self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.insert_kos += KEGG_list
                else:
                    self.no_insert_kos += KEGG_list
                # print (genome, gene_interval, KEGG_list, locate_insert_flag)
            # break
        print ("kos", len(self.insert_kos), len(self.no_insert_kos))
        return self.insert_kos, self.no_insert_kos

    def classify_bkp_kos(self, acc_file):
        bkp_dict = {}
        # self.near = surround


        my_bkps = self.read_bkp(acc_file)
        for bkp in my_bkps:
            if bkp.from_ref not in bkp_dict:
                bkp_dict[bkp.from_ref] = []
            bkp_dict[bkp.from_ref].append(bkp.from_bkp)
            if bkp.to_ref not in bkp_dict:
                bkp_dict[bkp.to_ref] = []
            bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "KEGG" not in gene_anno_dict:
                    continue
                if re.search("IS[0-9]", gene_anno_dict["product"]):  # IS element
                    continue
                KEGG_list = gene_anno_dict["KEGG"].split(",")
                new_KEGG_list = []
                for element in KEGG_list:
                    if element[:3] == "ko:":
                        element = element[3:]
                        new_KEGG_list.append(element)
                KEGG_list = new_KEGG_list

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    bkp_ko += KEGG_list
                else:
                    no_bkp_ko += KEGG_list
        # print_data(bkp_ko, bkp_ko_file)
        # print_data(no_bkp_ko, no_bkp_ko_file)
        return bkp_ko, no_bkp_ko
    
    def read_bkp(self, bkp_file):
        my_bkps = []
        f = open(bkp_file)
        all_rows = csv.reader(f)
        total_HGT_split_num = 0
        reads_num = 0
        for row in all_rows:
            if row[0][0] == "#":
                reads_num = int(row[0].split(";")[0].split(":")[1])
                # print (row, self.reads_num)
                pass
            elif row[0] == "from_ref":
                pass
            else:
                if reads_num == 0:
                    print ("old bkp", bkp_file)
                    # os.system("rm %s"%(bkp_file))
                    break
                eb = Acc_Bkp(row)
                eb.abundance = eb.cross_split_reads/reads_num
                if eb.from_ref_genome == eb.to_ref_genome:
                    continue
                if eb.abundance < abun_cutoff:
                    continue
                # if self.get_genome_taxa(eb.from_ref_genome) == self.get_genome_taxa(eb.to_ref_genome): # classify intra-genus and inter-genus HGTs
                #     continue
                total_HGT_split_num += eb.cross_split_reads
                # if eb.hgt_tag not in self.all_hgt:
                my_bkps.append(eb)
        for eb in my_bkps:
            eb.split_abundance = eb.cross_split_reads/total_HGT_split_num
        f.close()
        return my_bkps    

class Extract_COG(Extract_KO):

    def __init__(self, HGT_event_dict, surround=0):
        Extract_KO.__init__(self, HGT_event_dict, surround)
        self.transfer_cog = []
        self.no_transfer_cog = []
        self.insert_cog = []
        self.no_insert_cog = []
        self.bkp_cog = []
        self.no_bkp_cog = []
        self.data = []

    def classify_cog(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                # print (gene_anno_dict)
                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]
                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    for i in range(len(cog)):
                        self.transfer_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_transfer_cog += [cog[i]]
                # print (KEGG_list, locate_transfer_flag)
            # break
        print (len(self.transfer_cog), len(self.no_transfer_cog))
        # self.data = enrichment_analysis(self.transfer_cog, self.no_transfer_cog, "transfer", self.data)
        return self.transfer_cog, self.no_transfer_cog

    def classify_cog_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    for i in range(len(cog)):
                        self.insert_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_insert_cog += [cog[i]]

        print (len(self.insert_cog), len(self.no_insert_cog))
        # self.data = enrichment_analysis(self.insert_cog, self.no_insert_cog, "insert", self.data)
        return self.insert_cog, self.no_insert_cog

    def classify_bkp_cog(self, acc_file, surround):
        self.near = surround
        bkp_dict = {}


        my_bkps = self.read_bkp(acc_file)
        for bkp in my_bkps:
            if bkp.from_ref not in bkp_dict:
                bkp_dict[bkp.from_ref] = []
            bkp_dict[bkp.from_ref].append(bkp.from_bkp)
            if bkp.to_ref not in bkp_dict:
                bkp_dict[bkp.to_ref] = []
            bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]

                if "COG" not in gene_anno_dict:
                    continue
                cog = gene_anno_dict["COG"]

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    for i in range(len(cog)):
                        self.bkp_cog += [cog[i]]
                else:
                    for i in range(len(cog)):
                        self.no_bkp_cog += [cog[i]]
        # print (len(self.bkp_cog), len(self.no_bkp_cog))
        # self.data = enrichment_analysis(self.bkp_cog, self.no_bkp_cog, "BKP", self.data)
        return self.bkp_cog, self.no_bkp_cog

    def main(self):
        self.classify_regions()
        self.classify_bkp_cog()
        self.classify_cog()
        self.classify_cog_insert()

        df = pd.DataFrame(self.data, columns = ["category", "category_detail", "p_value", "fold", "gene_num", "profile", "locus_type"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        df.to_csv(cog_enrich, sep=',')
        print ("enriched COG num", len(self.data))

class Extract_product(Extract_KO):

    def __init__(self, HGT_event_dict):
        Extract_KO.__init__(self, HGT_event_dict)
        self.transfer_product = []
        self.no_transfer_product = []
        self.insert_product = []
        self.no_insert_product = []
        self.bkp_product = []
        self.no_bkp_product = []
        self.data = []

    def classify_product(self):  # gene in transfer region, or not
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                # print (gene_anno_dict)
                prod_category = annotation.classify_product(gene_anno_dict["product"])
                if prod_category == "unclassified":
                    continue

                #### check if the gene locates in transfer region
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                if locate_transfer_flag:
                    self.transfer_product += [prod_category]
                else:
                    self.no_transfer_product += [prod_category]
                # print (KEGG_list, locate_transfer_flag)
            # break
        print (len(self.transfer_product), len(self.no_transfer_product))
        self.data = enrichment_analysis_product(self.transfer_product, self.no_transfer_product, "transfer", self.data)

    def classify_product_insert(self):  # gene in insert site, or not
        for genome in self.insert_sites:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(self.insert_sites[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                prod_category = annotation.classify_product(gene_anno_dict["product"])

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.insert_product += [prod_category]
                else:
                    self.no_insert_product += [prod_category]

        print (len(self.insert_product), len(self.no_insert_product))
        self.data = enrichment_analysis_product(self.insert_product, self.no_insert_product, "insert", self.data)

    def classify_bkp_product(self):
        all_acc_file = hgt_result_dir + "/acc.list"
        os.system(f"ls {hgt_result_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        bkp_dict = {}

        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            for bkp in my_bkps:
                if bkp.from_ref not in bkp_dict:
                    bkp_dict[bkp.from_ref] = []
                bkp_dict[bkp.from_ref].append(bkp.from_bkp)
                if bkp.to_ref not in bkp_dict:
                    bkp_dict[bkp.to_ref] = []
                bkp_dict[bkp.to_ref].append(bkp.to_bkp)
        
        bkp_ko = []
        no_bkp_ko = []
        for genome in bkp_dict:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            insert_list = sorted(bkp_dict[genome])
            for gene_interval in gene_intervals:
                gene_anno_dict = annotation.gene_annotation[genome][str(gene_interval[0]) + "_" + str(gene_interval[1])]
                prod_category = annotation.classify_product(gene_anno_dict["product"])
                # if re.search("IS[0-9]", gene_anno_dict["product"]):  # IS element
                #     continue
                if remove_transposon_flag:
                    if annotation.classify_product(gene_anno_dict["product"]) == "Transposon":
                        continue

                #### check if the gene locates in insert site
                locate_insert_flag = False
                for site in insert_list:
                    if site > gene_interval[0] - self.near and site < gene_interval[1] + self.near:
                        locate_insert_flag = True
                        break
                if locate_insert_flag:
                    self.bkp_product += [prod_category]
                else:
                    self.no_bkp_product += [prod_category]
        print (len(self.bkp_product), len(self.no_bkp_product))
        self.data = enrichment_analysis_product(self.bkp_product, self.no_bkp_product, "bkp", self.data)

    def main(self):
        self.classify_regions()
        self.classify_bkp_product()
        self.classify_product()
        self.classify_product_insert()

        df = pd.DataFrame(self.data, columns = ["category", "p_value", "fold", "gene_num", "locus_type"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        df.to_csv(product_enrich, sep=',')
        # print ("enriched product num", len(self.data))

class Extract_seq(Extract_KO):

    def __init__(self, HGT_event_dict):
        Extract_KO.__init__(self, HGT_event_dict)
        self.ref_fasta = Fasta(database)
        # self.records_dict = SeqIO.index(database, "fasta")

    def transfer_cds(self):  # gene in transfer region, or not
        self.classify_regions()
        f = open(transfer_cds_fasta, 'w')
        h = open(no_transfer_cds_fasta, "w")

        cds_index = 0
        conserve_cds_index = 0
        for genome in self.transfer_regions:
            ### for each genome
            if genome not in annotation.gene_annotation:
                continue  # skip the genome without genes
            gene_intervals = annotation.gene_annotation[genome]["intervals"]
            transfer_intervals = self.transfer_regions[genome]
            for gene_interval in gene_intervals:
                
                locate_transfer_flag = check_overlap([gene_interval[0], gene_interval[1]], transfer_intervals, self.min_gene_frac)
                transfer_seq = self.ref_fasta[genome][gene_interval[0]:gene_interval[1]].seq
                # transfer_seq = self.records_dict[genome].seq[gene_interval[0]-1:gene_interval[1]]
                if genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) == "GUT_GENOME000269_2_220262_220672":
                    print ("Bug here", transfer_seq)
                if locate_transfer_flag:
                    cds_index += 1
                    cds_ID = genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) + "_" + str(cds_index) 
                    print (">%s"%(cds_ID), file = f)
                    print (transfer_seq, file = f)
                else:
                    conserve_cds_index += 1
                    cds_ID = genome + "_" + str(gene_interval[0]) + "_" + str(gene_interval[1]) + "_" + str(conserve_cds_index) 
                    print (">%s"%(cds_ID), file = h)
                    print (transfer_seq, file = h)                    
        f.close()
        h.close()

def enrichment_analysis(my_list, background_list, locus_type, data):
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
        print (category, p_value, oddsratio, a, b, c, d) 
        data.append([category, COG_dict[category], p_value, oddsratio, a, COG_profile_dict[category], locus_type])
    return data

def enrichment_analysis_product(my_list, background_list, locus_type, data):
    my_dict = Counter(my_list)
    background_dict = Counter(background_list)
    # data = []
    for prod_category in background_dict:
        if prod_category in my_dict:
            a = my_dict[prod_category]
        else:
            a = 0
        b = len(my_list) - a
        if prod_category in background_dict:
            c = background_dict[prod_category]
        else:
            c = 0
        d = len(background_list) - c

        oddsratio, p_value = fisher_exact([[a, b], [c, d]])
        print (locus_type, prod_category, p_value, oddsratio, a, b, c, d) 
        data.append([prod_category, p_value, oddsratio, a, locus_type])
    return data

def check_overlap(A, intervals, min_gene_frac):
    # Calculate the length of interval A
    A_length = A[1] - A[0]
    
    # Iterate through all the intervals in the list
    for interval in intervals:
        # Calculate the overlap between interval A and the current interval
        overlap = min(A[1], interval[1]) - max(A[0], interval[0])
        
        # If there is an overlap and the length of the overlap is longer than half of A, return True
        if overlap > 0 and overlap > A_length * min_gene_frac:
            return True
    
    # If there is no overlap that meets the criteria, return False
    return False

def print_data(my_set, file):
    f = open(file, 'w')
    for element in my_set:
        if element[:3] == "ko:":
            element = element[3:]
        print(element, end = "\n", file = f)
    f.close()   

class Classify(): # check the composition of transferred sequences

    def __init__(self, HGT_event_dict):
        self.HGT_event_dict = HGT_event_dict
        self.min_gene_frac = min_gene_frac

    def main(self):
        window = 100
        record_segment_tag = {}
        total_num, IS_num, trans_num = 0, 0, 0
        for sample in self.HGT_event_dict:
            for event in self.HGT_event_dict[sample]:
                segment_tag = "&".join([event.del_genome, str(round(event.del_start/window)), str(round(event.del_end/window))])  #, event.ins_genome_pure
                if segment_tag in record_segment_tag:
                    continue
                event.check_IS(self.min_gene_frac, annotation)
                if event.IS_flag:
                    IS_num += 1
                if event.Transposon_flag:
                    trans_num += 1
                total_num += 1
                record_segment_tag[segment_tag] = 1

            # break
        print (IS_num, total_num, IS_num/total_num, trans_num, trans_num/total_num, IS_num/trans_num)

def read_meta(meta_data):
    sra_sample_dict = {}
    for line in open(meta_data):
        if line.strip() == '':
            continue
        array = line.strip().split(',')
        if array[0] != 'Run':
            sra_id = array[0]
            sample_id = array[-2]
            if re.search("_", sample_id):
                sample_id = sample_id.split("_")[1]
            sra_sample_dict[sra_id] = sample_id
    return sra_sample_dict

def read_design():
    
    sample_individual_dict = {}
    sample_time_point = {}

    for line in open(design_file):
        if line.strip() == '':
            continue
        array = line.strip().split()
        if array[0] != 'Sample':
            sample_id = array[0]
            individual = array[3]
            sample_individual_dict[sample_id] = individual
            sample_time_point[sample_id] = int(array[4])
    return sample_individual_dict, sample_time_point

def get_individual_ID(): ## for the hybrid dataset
    sra2individual_ID = {}
    sra_sample_dict = read_meta(meta_data)
    sample_individual_dict, sample_time_point = read_design()
    for sra in sra_sample_dict:
        sample_id = sra_sample_dict[sra]
        individual_id = sample_individual_dict[sample_id]
        sra2individual_ID[sra] = individual_id
    return sra2individual_ID

def count_uniq_event_ratio(HGT_event_dict): # cal the ratio of uniq HGT events which occur only in one sample
    meta_data = "/mnt/d/HGT/time_lines/SRP366030.csv.txt"
    sra_sample_dict = read_meta(meta_data)

    uniq_event_dict = defaultdict(set)
    for sample in HGT_event_dict:
        for event in HGT_event_dict[sample]:
            if sample in sra_sample_dict:
                sample = sra_sample_dict[sample]
            uniq_event_dict[event.tag].add(sample)
    all_event_num = len(uniq_event_dict)
    uniq_event_num = 0
    for event_tag in uniq_event_dict:
        if len(uniq_event_dict[event_tag]) == 1:
            uniq_event_num += 1
    print (uniq_event_num, all_event_num, uniq_event_num/all_event_num)

def read_phenotype():
    phenotype_dict = {}
    pheno_result = "/mnt/d/HGT/association/phenotype.csv"
    for line in open(pheno_result):
        array = line.strip().split(",")
        ID = array[1]
        if ID == "sample":
            continue
        pheno = array[2:5]
        phenotype_dict[ID] = pheno
    return phenotype_dict


def load_gff_list(group_1_gffs):
    kos = []
    cogs = []
    for gff in group_1_gffs:
        genome_annotation = Annotation(gff) # load gene annotation
        kos += genome_annotation.kos
        # print ( genome_annotation.kos)
        cogs += genome_annotation.cogs

    for i in range(len(kos)):
        if kos[i][:3] == "ko:":
            kos[i] = kos[i][3:]

    return kos, cogs

def enrichment_run(group1_cogs, group2_cogs, group1_kos, group2_kos, cog_enrich, kegg_output):
    ## COG enrichment analysis
    data = enrichment_analysis(group1_cogs, group2_cogs, "cog", [])
    df = pd.DataFrame(data, columns = ["category", "category_detail", "p_value", "fold", "gene_num", "profile", "locus_type"])
    reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
    df["p.adj"] = pvals_corrected
    df.to_csv(cog_enrich, sep=',')
    print ("enriched COG num", len(data))


    ## pathway 
    print (len(group1_kos), group1_kos[:10])
    print (len(group2_kos))
    input_counts = kegg_enrichment.get_pathways(group1_kos, ko_pathway_dict)
    background_counts = kegg_enrichment.get_pathways(group2_kos, ko_pathway_dict)
    kegg_enrichment.enrichment_analysis(group1_kos, group2_kos, input_counts, background_counts, kegg_output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Functional analysis of HGT breakpoints", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--bkp", default = "NA", type=str, help="<str> concerned HGT breakpoints", metavar="\b")
    required.add_argument("--event", default = "NA", type=str, help="<str> concerned HGT events", metavar="\b")
    required.add_argument("--gff", type=str, help="<str> database gene annotation file", metavar="\b")
    # required.add_argument("--groupid", type=str, help="<str> column name used for grouping, default: phenotype", metavar="\b")
    # required.add_argument("--uhgg_meta", type=str, default = "genomes-all_metadata.tsv", help="<str> UHGG taxonomy metadata", metavar="\b")
    # required.add_argument("--kegg_output", type=str, default= "./pathway_enriched.csv", help="<str> KEGG output", metavar="\b")
    # required.add_argument("--cog_output", type=str, default= "./cog_enrich.csv",help="<str> COG output", metavar="\b")
    required.add_argument("--outdir", type=str, default= "./",help="<str> outdir", metavar="\b")
    required.add_argument("--out_prefix", type=str, default= "out",help="<str> out_prefix", metavar="\b")
    required.add_argument("--surround", type=int, default=5000, help="<str> focus on genes with distance within this value", metavar="\b")
    # required.add_argument("--group1", type=str, default = "D006262", help="<str> group1 name.", metavar="\b")
    # required.add_argument("--group2", type=str,  default = "D001249",  help="<str> group2 name.", metavar="\b")
    # required.add_argument("--breakpoint_dir", type=str, help="<str> HGT breakpoint_dir.", metavar="\b")
    # required.add_argument("-s", type=str, help="<str> split reads bam file.", metavar="\b")
    optional.add_argument("--abun_cutoff", type=float, default=1e-7, help="bkp abun cutoff", metavar="\b")
    optional.add_argument("--min_gene_frac", type=float, default=0.5, help="minimum gene fraction", metavar="\b")
    # optional.add_argument("-b", type=str, help="bed file of extracted ref.", metavar="\b")
    # optional.add_argument("-n", type=int, default=1, help="<0/1> 1 indicates the ref is extracted using kmer.", metavar="\b")
    # optional.add_argument("--read_info", type=int, default=1, help="<0/1> 1 indicates including reads info, 0 indicates not (just for evaluation).", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    abun_cutoff = args["abun_cutoff"]
    min_gene_frac = args["min_gene_frac"]
    gff = args["gff"]



    with open(sys.path[0] + '/ko_pathway_dict.pickle', 'rb') as f: ## this file is in meta_data/ko_pathway_dict.pickle
        ko_pathway_dict = pickle.load(f)
    # group_1_gffs = ["/mnt/d/breakpoints/HGT/UHGG/annotation/MGYG-HGUT-00823.gff", "/mnt/d/breakpoints/HGT/UHGG/annotation/MGYG-HGUT-02071.gff"]
    # group_2_gffs = ["/mnt/d/breakpoints/HGT/UHGG/annotation/MGYG-HGUT-00381.gff", "/mnt/d/breakpoints/HGT/UHGG/annotation/MGYG-HGUT-04193.gff"]
    COG_dict, COG_profile_dict = get_COG_dict()

    annotation = Annotation(gff) # load gene annotation
    annotation.read_gff()
    print ("annotation is done.")


    ### output
    # cog_enrich = args["cog_output"] # store COG enrichment results
    # kegg_output = args["kegg_output"] # to store enriched pathways
    out_prefix = args["out_prefix"]
    ## check if outdir exists, if not, create it
    if not os.path.exists(args["outdir"]):
        os.makedirs(args["outdir"])

    breakpoint_file = args["bkp"]

    if breakpoint_file != "NA":
        cog_enrich = args["outdir"] + "/" + out_prefix + "_COG_enrich.csv"
        kegg_output = args["outdir"] + "/" + out_prefix + "_KEGG_enrich.csv"
        ## only consider bkp
        extract = Extract_KO({}, args["surround"])
        group1_kos, group2_kos = extract.classify_bkp_kos(breakpoint_file)

        cog = Extract_COG({}, args["surround"])
        group1_cogs, group2_cogs = cog.classify_bkp_kos(breakpoint_file)
        enrichment_run(group1_cogs, group2_cogs, group1_kos, group2_kos, cog_enrich, kegg_output)

    identified_hgt = args["event"]
    if identified_hgt != "NA":


        cog_enrich = args["outdir"] + "/" + out_prefix + "_transfer_COG_enrich.csv"
        kegg_output = args["outdir"] + "/" + out_prefix + "_transfer_KEGG_enrich.csv"

        trans = Transfer_times()
        trans.read_events(identified_hgt)

        extract = Extract_KO(trans.HGT_event_dict, args["surround"])
        extract.classify_regions()
        group1_kos, group2_kos = extract.classify_kos()

        cog = Extract_COG(trans.HGT_event_dict, args["surround"])
        cog.classify_regions()
        group1_cogs, group2_cogs = cog.classify_cog()
        enrichment_run(group1_cogs, group2_cogs, group1_kos, group2_kos, cog_enrich, kegg_output)


        cog_enrich = args["outdir"] + "/" + out_prefix + "_insert_COG_enrich.csv"
        kegg_output = args["outdir"] + "/" + out_prefix + "_insert_KEGG_enrich.csv"

        trans = Transfer_times()
        trans.read_events(identified_hgt)

        extract = Extract_KO(trans.HGT_event_dict, args["surround"])
        extract.classify_regions()
        group1_kos, group2_kos = extract.classify_kos_insert()

        cog = Extract_COG(trans.HGT_event_dict, args["surround"])
        cog.classify_regions()
        group1_cogs, group2_cogs = cog.classify_cog_insert()
        enrichment_run(group1_cogs, group2_cogs, group1_kos, group2_kos, cog_enrich, kegg_output)






    








        
