#!/usr/bin/env python3

"""
Test the enrichment of microhomology in HGT breakpoint junctions
"""


import sys
import re, os
import csv
from scipy import stats
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
from pyfaidx import Fasta
from skbio.alignment import local_pairwise_align_ssw
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA
import re
from Bio.Seq import Seq
import pandas as pd
import random
import argparse

# from bkp_match import read_meta, read_design
# from bkp_match import Acc_Bkp

bin_size = 100
split_cutoff = 0  #10
sample_cutoff = 0  # 8
abun_cutoff = 1e-7  #1e-7

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
        self.score = float(list[11])

        self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))
        self.abundance = None
        self.split_abundance = None

        self.read = None

def read_samples(acc_file):
    all_bkps = []
    
    # for line in open(all_acc_file):
    #     acc_file = line.strip()
    
    my_bkps = read_bkp(acc_file)
    all_bkps += my_bkps

    print (len(all_bkps))
    return all_bkps
    
def read_bkp(bkp_file):
    my_bkps = []
    f = open(bkp_file)
    all_rows = csv.reader(f)
    sample_dict = {}
    total_HGT_split_num = 0
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
                break
            eb = Acc_Bkp(row)
            eb.abundance = eb.cross_split_reads/reads_num
            if eb.from_ref_genome == eb.to_ref_genome:
                continue
            if eb.cross_split_reads < split_cutoff:
                continue
            if eb.abundance < abun_cutoff:
                continue
            total_HGT_split_num += eb.cross_split_reads
            sample_dict[eb.hgt_tag] = 1
            my_bkps.append(eb)
    return my_bkps

def read_meta():
    
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


class Micro_homo():

    def __init__(self, workdir="/mnt/d/breakpoints/HGT/micro_homo/"):
        self.all_data, self.sample_num = None, None
        
        self.ref = database
        self.workdir = workdir
        self.ref_fasta = Fasta(self.ref)        
        self.cutoff = args["seq_len"]
        self.min_score = 8
        self.all_bkp = {}
        self.all_bkp_sample_num = {}
        self.tole_diff = 10
        self.shortest_len = 4
        # self.sam_number = 50000

        self.random_homo_seq_count = {}
        self.hgt_homo_seq_count = {}

    def get_reverse_complement_seq(self, sequence):
        sequence = sequence[::-1]
        trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
        string = sequence.translate(trantab)
        return string

    def extract_ref_seq(self, scaffold_name, start, end):
        if start < 1:
            start = 1
        # print (scaffold_name, start, end)
        return self.ref_fasta[scaffold_name][start:end].seq

    def ana_HGT_bkp(self):
        micro_freq_dict = {}
        homo_len_list = []
        index = 0
        hit_num = 0
        bkp_key_list = list(self.all_bkp.keys())
        random.shuffle(bkp_key_list)
        for bkp_key in bkp_key_list:
            bkp = self.all_bkp[bkp_key]

            if bkp.cross_split_reads < 5:
                continue

            max_homology_len = self.for_each_bkp(bkp)
            if max_homology_len == -1:
                continue
            if max_homology_len not in micro_freq_dict:
                micro_freq_dict[max_homology_len] = 0
            micro_freq_dict[max_homology_len] += 1
            homo_len_list.append(max_homology_len)
            index += 1
            if index % 1000 == 0:
                print (index, micro_freq_dict)  
            if index == self.sam_number:
                break          
        return micro_freq_dict, homo_len_list

    def ramdom_bkp_pair(self):
        if len(self.all_bkp) == 0:
            print ("no bkp?")  
        micro_freq_dict = {}
        homo_len_list = []
        index = 0
        hit_num = 0
        bkp_key_list = list(self.all_bkp.keys())
        random.shuffle(bkp_key_list)
        # for bkp_key in bkp_key_list:
    
        while True:
            i = np.random.randint(0, len(bkp_key_list))
            j = np.random.randint(0, len(bkp_key_list))
            if i == j:
                continue
            bkp1 = self.all_bkp[bkp_key_list[i]]
            bkp2 = self.all_bkp[bkp_key_list[j]]

            if bkp1.cross_split_reads < 5 or bkp2.cross_split_reads < 5:
                continue

            max_homology_len = self.for_two_random_bkps(bkp1, bkp2)
            if max_homology_len == -1:
                continue
            if max_homology_len not in micro_freq_dict:
                micro_freq_dict[max_homology_len] = 0
            micro_freq_dict[max_homology_len] += 1
            homo_len_list.append(max_homology_len)
            index += 1
            if index % 1000 == 0:
                print (index, micro_freq_dict)  
            if index == self.sam_number:
                break          
        return micro_freq_dict, homo_len_list

    def main(self):
        self.all_data, self.sample_num = get_data()
        print ("data is loaded")
        index = 0
        hit_num = 0
        sample_num = 0
        for sample in self.all_data:
            for bkp in sample.bkps:
                score = self.for_each_bkp(bkp, index)
                if score >= self.min_score:
                    hit_num += 1
                index += 1
                if index % 1000 == 0:
                    print (sample_num, index, hit_num, hit_num/index)
            sample_num += 1
            #     if index > 100:
            #         break
            # break
        print ("hit rate:", hit_num/index)

    def for_each_bkp(self, bkp):
        bkp.from_bkp -= 1
        from_seq = self.extract_ref_seq(bkp.from_ref, bkp.from_bkp-self.cutoff, bkp.from_bkp+self.cutoff)
        to_seq = self.extract_ref_seq(bkp.to_ref, bkp.to_bkp-self.cutoff, bkp.to_bkp+self.cutoff)
        if bkp.from_strand == "-":
            from_seq = self.get_reverse_complement_seq(from_seq)      
        if bkp.to_strand == "-":
            to_seq = self.get_reverse_complement_seq(to_seq)
        if re.search(">", from_seq) or re.search(">", to_seq):
            return -1, -1, '', ''
        if countN(from_seq) > 0 or countN(to_seq) > 0:
            return -1, -1, '', ''
        # test = ["GUT_GENOME018982_31", "GUT_GENOME239728_21"]
        # if bkp.from_ref in test and bkp.to_ref in test:
        #     print (bkp.from_ref, bkp.to_ref, from_seq, to_seq, self.find_mh(from_seq, to_seq))
        # score = self.get_mic_homo(from_seq, to_seq)
        # return score
        return self.get_micro_homo(from_seq, to_seq)

    def for_each_artifical_bkp(self, from_ref, from_bkp, to_ref, to_bkp, reverse_flag):
        from_bkp -= 1
        from_seq = self.extract_ref_seq(from_ref, from_bkp-self.cutoff, from_bkp+self.cutoff)
        to_seq = self.extract_ref_seq(to_ref, to_bkp-self.cutoff, to_bkp+self.cutoff)
        if bool(reverse_flag) == True:
            from_seq = self.get_reverse_complement_seq(from_seq)      

        if re.search(">", from_seq) or re.search(">", to_seq):
            return -1, -1, '', ''
        if countN(from_seq) > 0 or countN(to_seq) > 0:
            return -1, -1, '', ''
        # test = ["GUT_GENOME018982_31", "GUT_GENOME239728_21"]
        # if bkp.from_ref in test and bkp.to_ref in test:
        #     print (bkp.from_ref, bkp.to_ref, from_seq, to_seq, self.find_mh(from_seq, to_seq))
        # score = self.get_mic_homo(from_seq, to_seq)
        # return score
        return self.get_micro_homo(from_seq, to_seq)

    def for_two_random_bkps(self, bkp1, bkp2):
        bkp1.from_bkp -= 1
        from_seq = self.extract_ref_seq(bkp1.from_ref, bkp1.from_bkp - self.cutoff, bkp1.from_bkp + self.cutoff)
        if bkp1.from_strand == "-":
            from_seq = self.get_reverse_complement_seq(from_seq) 
        to_seq = self.extract_ref_seq(bkp2.to_ref, bkp2.to_bkp - self.cutoff, bkp2.to_bkp + self.cutoff)
        if bkp2.to_strand == "-":
            to_seq = self.get_reverse_complement_seq(to_seq)
        if re.search(">", from_seq) or re.search(">", to_seq):
            return -1
        if countN(from_seq) > 0 or countN(to_seq) > 0:
            return -1
        # test = ["GUT_GENOME018982_31", "GUT_GENOME239728_21"]
        # if bkp.from_ref in test and bkp.to_ref in test:
        #     print (bkp.from_ref, bkp.to_ref, from_seq, to_seq, self.find_mh(from_seq, to_seq))
        # score = self.get_mic_homo(from_seq, to_seq)
        # return score
        return self.get_micro_homo(from_seq, to_seq)

    def random_seq(self):
        micro_freq_dict = {}
        chroms = list(self.ref_fasta.keys())
        chroms_num = len(chroms)
        hit_num = 0
        # for i in range(1000):
        index = 0
        while True:
            seq_list = []
            for i in range(2):
                chrom = list(self.ref_fasta.keys())[np.random.randint(chroms_num)]
                chrom_len = len(self.ref_fasta[chrom])
                locus = np.random.randint(self.cutoff, chrom_len-self.cutoff)
                seq = self.extract_ref_seq(chrom, locus-self.cutoff, locus+self.cutoff)
                seq_list.append(seq)
                # print (chrom, locus)
            # score = self.get_mic_homo(seq_list[0], seq_list[1])
            if countN(seq_list[0]) > 0 or  countN(seq_list[1]) > 0:
                continue
            max_homology_len = self.get_micro_homo(seq_list[0], seq_list[1])
            if max_homology_len not in micro_freq_dict:
                micro_freq_dict[max_homology_len] = 0
            micro_freq_dict[max_homology_len] += 1
            index += 1
            if index % 1000 == 0:
                print (index, micro_freq_dict)
            if index == self.sam_number:
                break
        # print ("random hit rate", micro_freq_dict)
        return micro_freq_dict

    def visulaize(self, interval, match_seq):
        vis_seq = ''
        for i in range(0, interval[0]):
            vis_seq += "-"
        vis_seq += match_seq
        # for i in range(interval[0], interval[1]+1):
        #     vis_seq += "*"
        for i in range(interval[1]+1, self.cutoff*2):
            vis_seq += "-"
        return vis_seq

    def output_fasta(self, from_seq, to_seq, index):
        seq1 = SeqRecord(Seq(from_seq),
                        id="from_seq")
        seq2 = SeqRecord(Seq(to_seq),
                        id="to_seq")
        SeqIO.write(seq1, self.workdir+"from_seq_%s.fasta"%(index), "fasta")
        SeqIO.write(seq2, self.workdir + "to_seq_%s.fasta"%(index), "fasta")
        command = "blastn -evalue 0.05 -word_size 4 -gapopen 5 -gapextend 2 -reward 2 -penalty -3 -dust no\
            -query %s -subject %s -out %s"%(self.workdir+"from_seq_%s.fasta"%(index),\
             self.workdir + "to_seq_%s.fasta"%(index), self.workdir + "blast_%s.out"%(index))
        print (command)

    def find_mh(self, seq1, seq2, method):
        flag = False
        seq1 = seq1.lower()
        seq2 = seq2.lower()
        # print (self.shortest_len)
        for i in range(len(seq1)-self.shortest_len+1):
            search_seq = seq1[i:i+self.shortest_len]
            mat = re.search(search_seq, seq2)
            # print (search_seq, seq2)
            if mat:
                if method == "random":
                    if search_seq not in self.random_homo_seq_count:
                        self.random_homo_seq_count[search_seq] = 0
                    self.random_homo_seq_count[search_seq] += 1

                if method == "hgt":
                    if search_seq not in self.hgt_homo_seq_count:
                        self.hgt_homo_seq_count[search_seq] = 0
                    self.hgt_homo_seq_count[search_seq] += 1

                mat_locus = mat.span()[0]
                # print (abs(mat_locus - i), self.tole_diff, "true")
                if abs(mat_locus - i) <= self.tole_diff :
                    # print (abs(mat_locus - i), "true")
                    flag =  True
        return flag

    def get_micro_homo(self, seq1, seq2):
        from_seq = DNA(seq1)
        to_seq = DNA(seq2)
        alignment, score, start_end_positions = global_pairwise_align_nucleotide(from_seq, to_seq) #, mismatch_score = -1000
        # alignment, score, start_end_positions = local_pairwise_align_ssw(from_seq, to_seq)

        max_homology_len = extract_homology(alignment)
        # if max_homology_len > 0:
        #     print ("fr", alignment[0])
        #     print ("to", alignment[1])
        return max_homology_len, alignment, from_seq, to_seq

def cal_ave_homo_len(homo_freq):
    sample_num = sum(list(homo_freq.values()))
    average_homo_len = 0
    for length in  homo_freq:
        average_homo_len += length * (homo_freq[length]/sample_num)
    return average_homo_len

def microhomology_freq_compare():
    data = []
    # for diff in range(1):
    #     for score in range(3, 11):
    mic = Micro_homo()
    mic.sam_number = 2000
    # mic.tole_diff = diff
    # mic.shortest_len = score
    # print ("score is", score, "diff is", diff)
    hgt_freq = mic.ana_HGT_bkp()
    print ("HGT", hgt_freq, cal_ave_homo_len(hgt_freq))
    # random_freq = mic.random_seq()
    # print ("random", random_freq)
    random_hgt_freq = mic.ramdom_bkp_pair()
    print ("random HGT", random_hgt_freq, cal_ave_homo_len(random_hgt_freq))
    # data.append([diff, score, hgt_freq, "HGT"])
    # data.append([diff, score, random_freq, "Random"])
    print ("-----------------------")
    # df = pd.DataFrame(data, columns = ["diff", "length", "frequency", "group"])
    # df.to_csv('/mnt/d/breakpoints/script/analysis/microhomo_freq.csv', sep='\t')

def count_seq_freq(mic, hgt_homo_seq_count):
    new_dict = {}
    print (len(hgt_homo_seq_count))
    for key in sorted(hgt_homo_seq_count):
        # print (key)
        if key not in new_dict:
            if mic.get_reverse_complement_seq(key) in new_dict:
                new_dict[mic.get_reverse_complement_seq(key)] += hgt_homo_seq_count[key]
            else:
                new_dict[key] = hgt_homo_seq_count[key]
        else:
            new_dict[key] += hgt_homo_seq_count[key]
            # print ("impossible")
    print (len(new_dict))
    return new_dict
        
def countN(sequence):
    # initialize a counter variable
    count = 0

    # loop through the sequence and count the number of "N" characters
    for char in sequence:
        if char.upper() == 'N':
            count += 1
    return count

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

def extract_homology(alignment):
    homology_list = []
    
    homology_flag = False
    homology_seq = ''
    for i in range(len(alignment[0])):
        if str(alignment[0][i]) != "-" and str(alignment[1][i]) != "-":
            # print (i, alignment[0][i], alignment[1][i])
            if homology_flag == False:
                homology_flag = True
                homology_seq += str(alignment[0][i])
            else:
                homology_seq += str(alignment[0][i])
        else:
            if len(homology_seq) > 0:
                homology_list.append(homology_seq)
            homology_flag = False
            homology_seq = ''
    if len(homology_seq) > 0:
        homology_list.append(homology_seq)  
    homology_len_list = [len(x) for x in homology_list]
    # if len(homology_len_list) > 0 and max(homology_len_list) > 20:
    #     print (homology_list)
    #     print (alignment)
    # print (homology_list)
    # return homology_list
    if len(homology_len_list) > 0:
        return max(homology_len_list)
    else:
        return 0

def length_cutoff_count(raw_len_list, cutoff):
    new_list = []
    for length in raw_len_list:
        if length > cutoff:
            new_list.append(1)
        else:
            new_list.append(0)
    return new_list, sum(new_list)/len(new_list)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Functional analysis of HGT breakpoints", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--bkp", default = "NA", type=str, help="<str> concerned HGT breakpoints", metavar="\b")
    required.add_argument("--event", default = "NA", type=str, help="<str> concerned HGT events", metavar="\b")
    required.add_argument("--database", type=str, help="<str> UHGG database file", metavar="\b")
    required.add_argument("--min_len", type=int,default=2, help="minimum microhomology length", metavar="\b")
    # required.add_argument("--uhgg_meta", type=str, default = "genomes-all_metadata.tsv", help="<str> UHGG taxonomy metadata", metavar="\b")
    required.add_argument("--output", type=str, default= "./micro.csv", help="<str> output", metavar="\b")
    # required.add_argument("--cog_output", type=str, default= "./cog_enrich.csv",help="<str> COG output", metavar="\b")
    # required.add_argument("--surround", type=int, default=5000, help="<str> focus on genes with distance within this value", metavar="\b")
    optional.add_argument("--seq_len", type=int, default=20, help="upstream and downstream length from the BKP pos", metavar="\b")
    # optional.add_argument("--min_gene_frac", type=float, default=0.5, help="minimum gene fraction", metavar="\b")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    
    # database = "/mnt/d/breakpoints/HGT/micro_homo/UHGG_reference.formate.fna" # UHGG v1 reference

    # load data
    database = args["database"]

    if args["bkp"] != "NA":
        all_bkps = read_samples(args["bkp"])
        # abun_cutoff = args["abun_cutoff"]
        # bin_size = args["bin_size"]

        data = []
        mic = Micro_homo()
        for bkp in all_bkps:
            
            max_homology_len, alignment, from_seq, to_seq = mic.for_each_bkp(bkp)
            # if max_homology_len == -1:
            #     continue
            # if max_homology_len < args["min_len"]:
            #     continue
            # print (bkp.from_ref)
            print (max_homology_len, alignment[0],  alignment[1], from_seq, to_seq, sep="\n")

            data.append(["sample", bkp.from_ref, bkp.from_bkp, bkp.from_strand, bkp.from_bkp-args["seq_len"], bkp.from_bkp+args["seq_len"],from_seq, \
                bkp.to_ref, bkp.to_bkp, bkp.to_strand, bkp.to_bkp-args["seq_len"], bkp.to_bkp+args["seq_len"],to_seq, max_homology_len, alignment[0],  alignment[1]])
        
        df = pd.DataFrame(data, columns = ["sample", "from_ref", "from_bkp", "from_strand", "from_start", "from_end", "from_seq",  \
        "to_ref", "to_bkp", "to_strand", "to_start", "to_end", "to_seq", "homology_len", "from_alignment", "to_alignment"])
        df.to_csv(args["output"], sep=',')  


    identified_hgt = args["event"]
    if identified_hgt != "NA":
        trans = Transfer_times()
        trans.read_events(identified_hgt)
        data = []
        mic = Micro_homo()
        for sample in trans.HGT_event_dict:
            for event in trans.HGT_event_dict[sample]:
                if event.reverse_flag == False:
                    from_ref = event.ins_genome
                    from_bkp = event.ins_pos
                    to_ref = event.del_genome
                    to_bkp = event.del_start
                    max_homology_len, alignment, from_seq, to_seq = mic.for_each_artifical_bkp(from_ref, from_bkp, to_ref, to_bkp, event.reverse_flag)
                    data.append([sample, from_ref, from_bkp, event.reverse_flag, from_bkp-args["seq_len"], from_bkp+args["seq_len"],from_seq, \
                to_ref, to_bkp, event.reverse_flag, to_bkp-args["seq_len"], to_bkp+args["seq_len"],to_seq, max_homology_len, alignment[0],  alignment[1]])

                    from_ref = event.del_genome
                    from_bkp = event.del_end
                    to_ref = event.ins_genome
                    to_bkp = event.ins_pos
                    max_homology_len, alignment, from_seq, to_seq = mic.for_each_artifical_bkp(from_ref, from_bkp, to_ref, to_bkp, event.reverse_flag)
                    data.append([sample, from_ref, from_bkp, event.reverse_flag, from_bkp-args["seq_len"], from_bkp+args["seq_len"],from_seq, \
                to_ref, to_bkp, event.reverse_flag, to_bkp-args["seq_len"], to_bkp+args["seq_len"],to_seq, max_homology_len, alignment[0],  alignment[1]])
                else:
                    from_ref = event.ins_genome
                    from_bkp = event.ins_pos
                    to_ref = event.del_genome
                    to_bkp = event.del_end
                    max_homology_len, alignment, from_seq, to_seq = mic.for_each_artifical_bkp(from_ref, from_bkp, to_ref, to_bkp, event.reverse_flag)
                    data.append([sample, from_ref, from_bkp, event.reverse_flag, from_bkp-args["seq_len"], from_bkp+args["seq_len"],from_seq, \
                to_ref, to_bkp, event.reverse_flag, to_bkp-args["seq_len"], to_bkp+args["seq_len"],to_seq, max_homology_len, alignment[0],  alignment[1]])
                    from_ref = event.del_genome
                    from_bkp = event.del_start
                    to_ref = event.ins_genome
                    to_bkp = event.ins_pos
                    max_homology_len, alignment, from_seq, to_seq = mic.for_each_artifical_bkp(from_ref, from_bkp, to_ref, to_bkp, event.reverse_flag)
                    data.append([sample, from_ref, from_bkp, event.reverse_flag, from_bkp-args["seq_len"], from_bkp+args["seq_len"],from_seq, \
                to_ref, to_bkp, event.reverse_flag, to_bkp-args["seq_len"], to_bkp+args["seq_len"],to_seq, max_homology_len, alignment[0],  alignment[1]])
                    # print (max_homology_len, alignment[0],  alignment[1], from_seq, to_seq, sep="\n")

        df = pd.DataFrame(data, columns = ["sample", "from_ref", "from_bkp", "from_strand", "from_start", "from_end", "from_seq",  \
        "to_ref", "to_bkp", "to_strand", "to_start", "to_end", "to_seq", "homology_len", "from_alignment", "to_alignment"])
        df.to_csv(args["output"], sep=',')  






    # compare the microhomology length between HGT breakpoint pairs and random two HGT breakpoints 
    # mic = Micro_homo()
    # mic.sam_number = 10000
    # all_bkps_dict = {}
    # for bkp in all_bkps:
    #     all_bkps_dict[bkp.hgt_tag] = bkp
    # mic.all_bkp = all_bkps_dict
    # hgt_freq, hgt_homo = mic.ana_HGT_bkp()
    # print ("HGT", sorted(hgt_freq), cal_ave_homo_len(hgt_freq))
    # random_hgt_freq, random_homo = mic.ramdom_bkp_pair()
    # U1, p = mannwhitneyu(hgt_homo, random_homo)
    # print ("random HGT", sorted(random_hgt_freq), cal_ave_homo_len(random_hgt_freq))
    # print (U1, p)
    # print ("-----------------------")

    # data = []
    # max_len = max(list(hgt_freq.keys()))
    # for length in range(max_len + 1):
    #     if length in hgt_freq:
    #         data.append([length, hgt_freq[length]/mic.sam_number, "HGT"])
    #     if length in random_hgt_freq:
    #         data.append([length, random_hgt_freq[length]/mic.sam_number, "Random"])   
    # df = pd.DataFrame(data, columns = ["length", "frequency", "group"])
    # df.to_csv('/mnt/d/breakpoints/script/analysis/microhomo_freq.csv', sep='\t')  

    # data = []
    # for i in range(mic.sam_number):
    #     data.append([hgt_homo[i], "HGT_Junction"])
    #     data.append([random_homo[i], "Expected"])
    # df = pd.DataFrame(data, columns = ["Microhomology", "Group"])
    # df.to_csv('/mnt/d/R_script_files//microhomology_length.csv', sep=',')  # the file to make figure

    # cutoff = 5
    # cut_hgt_homo, hgt_ratio = length_cutoff_count(hgt_homo, cutoff)
    # cut_random_homo, random_ratio = length_cutoff_count(random_homo, cutoff)
    # U1, p = mannwhitneyu(cut_hgt_homo, cut_random_homo)
    # print (U1, p, hgt_ratio, random_ratio)
    # print ("-----------------------")


      