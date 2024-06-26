from util import *
import sys
import subprocess
import getopt
import os
import argparse
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
import re, os
import csv

level_list = ["phylum", "class", "order", "family", "genus", "species", "genome"]

class Taxonomy():
    def __init__(self, UHGG_meta):
        self.taxonomy_dict = {}
        self.read_UHGG(UHGG_meta)
        
    def read_UHGG(self, UHGG_meta):
        df = pd.read_table(UHGG_meta) 
        for index, row in df.iterrows():
            # self.read_UHGG[row["Genome"]] = row["Species_rep"]
            genome = row["Genome"]
            lineage = row["Lineage"]
            self.taxonomy_dict[genome] = lineage

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
        self.abundance = None
        self.split_abundance = None
        self.read = None

        self.from_ref_lineage = taxonomy.taxonomy_dict[self.from_ref_genome]
        self.to_ref_lineage = taxonomy.taxonomy_dict[self.to_ref_genome]

        taxa1 = self.from_ref_lineage.split(";")[level]
        taxa2 = self.to_ref_lineage.split(";")[level]
        if taxa1[1:] == "__" or taxa2[1:] == "__":
            self.hgt_tag = "NA"
        else:
            self.hgt_tag = "&".join(sorted([taxa1, taxa2]))
        # self.hgt_tag = self.from_ref + "&" + str(int(self.from_bkp/bin_size)) + "&" + self.to_ref + "&" + str(int(self.to_bkp/bin_size))

class Sample():

    def __init__(self, bkp_list, ID, pheno):
        self.bkps = bkp_list
        self.ID = ID  
        # self.cohort = pheno[0]  
        self.disease = pheno 
        # self.full_disease = pheno[2].split(";")


class Data_load():

    def __init__(self, breakpoint_dir):
        self.sample_obj_list = []
        self.breakpoint_dir = breakpoint_dir
        
    def read_samples(self):

        all_acc_file = self.breakpoint_dir + "/acc.list"
        os.system(f"ls {self.breakpoint_dir}/*acc.csv |grep -v repeat >{all_acc_file}")
        # os.system(f"ls {tgs_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        # os.system(f"ls {wenkui_dir}/*acc.csv |grep -v repeat >>{all_acc_file}")
        
        for line in open(all_acc_file):
            acc_file = line.strip()
            sra_id = acc_file.split("/")[-1].split(".")[0]
            my_bkps = self.read_bkp(acc_file)
            print (sra_id, len(my_bkps))
            if len(my_bkps) > 0 and sra_id in phenotype_dict:
                sample = Sample(my_bkps, sra_id, phenotype_dict[sra_id])
                self.sample_obj_list.append(sample)
        print ("data is loaded.", "sample num:", len(self.sample_obj_list))
        
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


def get_genome_taxa(genome, level):
    # g1 = get_pure_genome(genome)
    taxa = taxonomy.taxonomy_dict[genome].split(";")[level]
    return taxa

def get_tag(bkp, level):
    if level == 7:
        from_tax = bkp.from_ref_genome
        to_tax = bkp.to_ref_genome
    else:
        from_tax = bkp.from_ref_lineage.split(";")[level] #  bkp.from_ref
        to_tax = bkp.to_ref_lineage.split(";")[level]  #bkp.to_ref
    from_tax = "_".join(from_tax.split())
    to_tax = "_".join(to_tax.split())
    tax = sorted([from_tax, to_tax])
    new_tag = "&".join(tax)       
    node1 = from_tax  
    node2 = to_tax  
    return new_tag


class Marker():

    def __init__(self, group1, group2, sample_obj_list):
        self.group1 = group1
        self.group2 = group2
        self.all_HGTs = None
        self.sample_count = None 
        self.marker_num = 100000000000
        self.markers = None
        self.split_data_dict = None
        self.sample_obj_list = sample_obj_list

    def extract_HGT(self, selected_samples):
        self.all_HGTs = {}
        self.sample_count = {self.group1:0, self.group2:0}
        sample_num = 0
        for sample in self.sample_obj_list:
            if sample.ID not in selected_samples:# self.split_data_dict:
                continue
            # elif self.split_data_dict[sample.ID] == "validation": #  Feature ranking was performed internally to each training fold to avoid overfitting.
            #     continue
            sample_num += 1
            self.sample_count[sample.disease] += 1
            sample_dict = {}
            for bkp in sample.bkps:
                if bkp.hgt_tag == "NA":
                    continue
                if bkp.hgt_tag in sample_dict:
                    continue
                if bkp.hgt_tag not in self.all_HGTs:
                    self.all_HGTs[bkp.hgt_tag] = {self.group1:0, self.group2:0}
                sample_dict[bkp.hgt_tag] = 1
                self.all_HGTs[bkp.hgt_tag][sample.disease] += 1
        
        # print (sample_num, self.sample_count)

        filtered_HGT = {}
        # print ("Bkp num in the two groups", len(self.all_HGTs))
        for hgt_tag in self.all_HGTs:
            # if (self.all_HGTs[hgt_tag][self.group1] + self.all_HGTs[hgt_tag][self.group2])/(self.sample_count[self.group1]+self.sample_count[self.group2]) < cutoff:
            #     pass
            # else:
            filtered_HGT[hgt_tag] = self.all_HGTs[hgt_tag]
        self.all_HGTs = filtered_HGT
        # print ("Filtered bkp num in the two groups", len(self.all_HGTs))
        # print ("%s num: %s, %s num %s."%(self.group1, self.sample_count[self.group1],self.group2, self.sample_count[self.group2]))

    def select_diff_HGT(self, output):
        hgt_p_value_dict = {}
        data = []
        for hgt_tag in self.all_HGTs:

            a = self.all_HGTs[hgt_tag][self.group1]
            b = self.sample_count[self.group1] - self.all_HGTs[hgt_tag][self.group1]
            c = self.all_HGTs[hgt_tag][self.group2]
            d = self.sample_count[self.group2] - self.all_HGTs[hgt_tag][self.group2]
            # print (a, b, c, d)
            # group1_freq = a/(a+b)
            # group2_freq = c/(c+d)
            # print (self.group1, self.sample_count[self.group1], self.group2, self.sample_count[self.group2], a, b, c, d)
            oddsratio, p_value = fisher_exact([[a, b], [c, d]])
            data.append([hgt_tag, p_value, oddsratio, a])

        df = pd.DataFrame(data, columns = ["genus_pair", "p_value", "oddsratio", "gp_num"])
        reject, pvals_corrected, _, alphacBonf = multipletests(list(df["p_value"]), alpha=0.05, method='bonferroni')
        df["p.adj"] = pvals_corrected
        df.to_csv(output, sep=',')

        # filtered_df = df[df['p.adj'] < 0.05]
        # filtered_df.to_csv(output, sep='\t')
        # for index, row in filtered_df.iterrows():
        #     hgt_p_value_dict[row["genus_pair"]] = row["p.adj"]

        # sorted_hgt_p_value = sorted(hgt_p_value_dict.items(), key=lambda item: item[1], reverse = False)
        # self.markers = {}
        # if len(hgt_p_value_dict) < self.marker_num:
        #     real_marker_num = len(hgt_p_value_dict)
        #     print ("***only has %s features."%(real_marker_num))
        # else:
        #     real_marker_num = self.marker_num

        # for i in range(real_marker_num):
        #     marker = sorted_hgt_p_value[i][0]
        #     p = sorted_hgt_p_value[i][1]
        #     self.markers[marker] = i
        #     # print ("marker", marker, p)
        # print (f"total marker number is {len(filtered_df)}, used biomarkers {real_marker_num}")
        # return len(self.markers) # total marker number

    def select_sample(self):
        selected_samples = []
        group1_num, group2_num = 0, 0
        for sample in self.sample_obj_list:

            # if len(sample.full_disease) != 1:
            #     continue
            # if sample.disease == "control" and sample.full_disease[0] != "healthy":
            #     continue
            if sample.disease == '':
                continue
            if sample.disease ==  self.group1:
                index = 0
                group1_num += 1
            elif sample.disease ==self.group2:
                index = 1
                group2_num += 1
            else:
                print ("discard group")
                continue 
            selected_samples.append(sample.ID)
        print (self.group1, group1_num, self.group2, group2_num)
        # shuffle(selected_samples)    
        return selected_samples

    
        print ("biomarker num:", len(X[0]))
        n_splits = 5
        cv = StratifiedKFold(n_splits=n_splits)
        # classifier = svm.SVC(kernel="linear", probability=True, random_state=self.random_state)
        classifier = RandomForestClassifier(n_estimators=100)

        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        fig, ax = plt.subplots(figsize=(6, 6))
        for fold, (train, test) in enumerate(cv.split(X, y)):
            classifier.fit(X[train], y[train])
            viz = RocCurveDisplay.from_estimator(
                classifier,
                X[test],
                y[test],
                name=f"ROC fold {fold}",
                alpha=0.3,
                lw=1,
                ax=ax,
                plot_chance_level=(fold == n_splits - 1),
            )
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)

        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        ax.plot(
            mean_fpr,
            mean_tpr,
            color="b",
            label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
            lw=2,
            alpha=0.8,
        )

        std_tpr = np.std(tprs, axis=0)
        tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
        tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
        ax.fill_between(
            mean_fpr,
            tprs_lower,
            tprs_upper,
            color="grey",
            alpha=0.2,
            label=r"$\pm$ 1 std. dev.",
        )

        ax.set(
            xlim=[-0.05, 1.05],
            ylim=[-0.05, 1.05],
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=f"Mean ROC curve with variability\n({self.group1} vs. {self.group2}, No. of biomarkers: {len(X[0])})",
        )
        ax.axis("square")
        ax.legend(loc="lower right")
        # plt.show()
        plt.savefig(f'/mnt/d/HGT/biomarker/class_{self.group1}_vs_{self.group2}.pdf')

def split_input_acc(all_acc, output_csv_path):
    ## check if output_csv_path exists, if not create it
    if not os.path.exists(output_csv_path):
        os.makedirs(output_csv_path)
    # Open the original CSV file
    with open(all_acc, mode='r', newline='') as file:
        for line in file:
            if re.search("#sample", line):
                sample_name = line.strip().split(":")[1]
                output_csv = f'{output_csv_path}/{sample_name}.acc.csv'
                out_f = open(output_csv, "w")
            else:
                print (line, file=out_f, end = '')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Identify HGT biomarkers between two groups", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--ann", type=str, help="<str> input file of group info", metavar="\b")
    required.add_argument("--groupid", type=str, help="<str> column name used for grouping, default: phenotype", metavar="\b")
    required.add_argument("--uhgg_meta", type=str, help="<str> UHGG taxonomy metadata, genomes-all_metadata.tsv", metavar="\b")
    # required.add_argument("--uhgg_meta", type=str, default = "genomes-all_metadata.tsv", help="<str> UHGG taxonomy metadata", metavar="\b")
    required.add_argument("--output", type=str, default = "HGT_biomarker.tsv", help="<str> HGT biomarker output", metavar="\b")
    required.add_argument("--taxa_level", type=int, default=5, help="<str> Original Metagenomic reference", metavar="\b")
    required.add_argument("--group1", type=str, default = "D006262", help="<str> group1 name.", metavar="\b")
    required.add_argument("--group2", type=str,  default = "D001249",  help="<str> group2 name.", metavar="\b")
    required.add_argument("--all_acc", type=str, help="<str> a file contains breakpoints of all samples.", metavar="\b")
    required.add_argument("--breakpoint_dir", type=str, help="<str> HGT breakpoint_dir to store the acc.csv file for each sample.", metavar="\b")
    # required.add_argument("-s", type=str, help="<str> split reads bam file.", metavar="\b")
    optional.add_argument("--abun_cutoff", type=float, default=1e-7, help="bkp abun cutoff", metavar="\b")
    # optional.add_argument("-b", type=str, help="bed file of extracted ref.", metavar="\b")
    # optional.add_argument("-n", type=int, default=1, help="<0/1> 1 indicates the ref is extracted using kmer.", metavar="\b")
    # optional.add_argument("--read_info", type=int, default=1, help="<0/1> 1 indicates including reads info, 0 indicates not (just for evaluation).", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    metadata = args["ann"]
    groupid = args["groupid"]
    level = args["taxa_level"]
    group1 = args["group1"]
    group2 = args["group2"]
    abun_cutoff = args["abun_cutoff"]
    output = args["output"]

    group = metadata2gf(metadata,groupid)

    # group.set_axis(group.iloc[0], axis=1, inplace=True)
    # group = group[1:] 

    phenotype_dict = {}
    for index, row in group.iterrows():
        # print (index, row.values[0])
        phenotype_dict[index] = row.values[0]

    taxonomy = Taxonomy(args["uhgg_meta"])
    # print (group)

    split_input_acc(args["all_acc"], args["breakpoint_dir"])
    dat = Data_load(args["breakpoint_dir"])
    dat.read_samples()

    biomarker = Marker(group1, group2, dat.sample_obj_list)
    selected_samples = biomarker.select_sample()
    biomarker.extract_HGT(selected_samples)
    biomarker.select_diff_HGT(output)



# '''
#     This is alpha_diversity.
#     options:
#     --abdf   <str> input file path
#     --dsf    <str> input file from subplatform
#     --ann    <str> input file of group info
#     --groupid  <str> column name used for grouping, default: phenotype
#     --method <str> method for alpha diversity, default: shannon. [shannon/simpson/invsimpson/ACE/Chao1/observedSpecies]
#     --reads  <int> number of reads used to calculate richness, need with methods: ACE/Chao1/observedSpecies, default: 500000
#     --testing <str> method for testing, default: t.test. [wilcox.test/t.test/kruskal.test/aov]
# '''

# python_file = os.path.abspath(__file__) 
# python_dir = os.path.dirname(python_file)
# script_path = os.path.join(python_dir, 'alpha.R')

# ifile1 = ''
# ifile2 = ''

# ops, args = getopt.getopt(sys.argv[1:], '', ['abdf=', 'dsf=', 'ann=', 'method=', 'reads=', 'groupid=',"testing="])
# for op, arg in ops:
#     if op == '--abdf':
#         ifile1 = arg
#     if op == '--dsf':
#         ifile2 = arg
#     if op == '--ann':
#         metadata = arg
#     if op == '--groupid':
#         groupid = arg
#     if op == '--method':
#         method = arg
#     if op == '--reads':
#         reads = arg
#     if op == '--testing':
#         testing = arg

# # merged_df = get_merged(ifile1,ifile2)
# group = metadata2gf(metadata,groupid)
# # if not check_valid(group, merged_df):
# #     exit(2)
# group_file = 'merged_input.group_info.tsv'
# group.to_csv(group_file, sep='\t')


# merged_df.to_csv('merged_input.abundance.all.tsv', sep='\t')
# split_tax = tax_split(merged_df)
# tax_list = [i for i in split_tax.keys()]

# output1 = 'output.alpha_diversity.tsv'
# output2 = 'output.pvalue.tsv'

# for i in range(len(tax_list)):
#     tax = tax_list[i]
#     split_tax[tax].index.name = tax
#     tax_abd = 'merged_input.abundance.' + tax +'.tsv'
#     split_tax[tax].to_csv(tax_abd, sep='\t')
#     if i == 0:
#         result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-r',reads,'-e',testing,'-a','F'],stdout=subprocess.PIPE)
#     else:
#         result = subprocess.run(['Rscript',script_path,'-i',tax_abd,'-g',group_file,'-o',output1,'-p',output2,'-t',tax,'-m',method,'-r',reads,'-e',testing,'-a','T'],stdout=subprocess.PIPE)
#     if result.returncode > 0:
#         print (result.stderr)
#         exit(1)
