import pandas as pd 
import numpy as np
import math
from shutil import copyfile
from scipy.stats import ttest_ind
import chardet
import numpy as np
from sys import exit
from collections import defaultdict





# version 1.0 author: Lijia Che
# initial implementation

# version 1.1 author: Yiqi Jiang
# add warning log, genarated group information from metadata, convert mpa to newick

# version 1.2 author: Yiqi Jiang
# add filter_metadata_att

    
def delete_zero_row(data):    
    # print filtered taxonomy
    filtered_data = data[data.apply(np.sum,axis = 1) == 0]
    for i in filtered_data._stat_axis.values.tolist():
        print('Warning: Discard taxonomy {} for no sample have profile in the profiling table.'.format(i))    
    
    data = data[data.apply(np.sum,axis = 1) != 0]
    tax = data._stat_axis.values.tolist()
    tax.sort()
    data = data.reindex(tax)
    return data

def tax_split(data1):
    samples = data1.columns.values.tolist()
    abd_type = data1._stat_axis.values.tolist()

    # kpcofgst 8 layers
    len_relation = [None, 'k', 'p', 'c', 'o', 'f', 'g', 's', 't']
    row_names = {
        'k': [], 'p': [], 'c': [], 'o': [], 'f':[], 'g':[], 's': [], 't': []
    }

    for t in abd_type:
        layers = t.split('|')
        nlayers = len(layers)
        row_names[len_relation[nlayers]].append(t)

    # split
    merged_data = {}
    for i in range(1,9):
        layer = len_relation[i]
        row_name = row_names[layer]
        data = pd.DataFrame(index=row_name, columns=samples)
        for type_name in row_name:
            for sample in samples:
                data.loc[type_name, sample] = data1.loc[type_name, sample]   
        merged_data[layer] = data.copy()
    if merged_data['t'].empty:
        del merged_data['t']
    return merged_data



# to check input 
def check_valid(gf, mergef):
    abd_slist = list(mergef.columns)
    g_slist = list(gf.index)
    valid = True
    for s in g_slist:
        if s not in abd_slist:
            print('Error: profile of sample {} in group information dose NOT exist.'.format(s))
            valid = False
    for s in list(abd_slist):
        if s not in g_slist:
            print('Error: group information of sample {} in profile dose NOT exist.'.format(s))
            valid = False
            
    # check all-0 sample
    if valid:
        allzero_sample = mergef.loc[:, (mergef == 0).all(axis=0)]
        if not allzero_sample.empty:
            for i in allzero_sample._stat_axis.values.tolist():
                print('Error: the profiling sum of sample {} is 0.'.format(i))
                valid = False
        
    return valid

def auto_decode(ifile):
    ec_type = chardet.detect(open(ifile, 'rb').read())['encoding']
    data = pd.read_csv(ifile, header=0, sep='\t', index_col = 0,encoding=ec_type)
    return data


    
def metadata2gf(metadata,group_id='phenotype'):
    data = auto_decode(metadata)
    if group_id in data.columns:
        gf = data[[group_id]]
        if len(data[group_id].unique()) != 2:
            print('Error: the column {} does not have exact 2 level.'.format(group_id))
            exit(1)
        return gf
    else:
        print('Error: the column {} does not exists in the metadata input.'.format(group_id))
        exit(1)

# main merge function
def get_merged(ifile1,ifile2):
    if ifile1 in ['', None, 'None'] and ifile2 in ['', None, 'None']:
        print('Error: please input profile.')
        exit(1)
    elif ifile1 in ['', None, 'None']:
        merged_df = auto_decode(ifile2)
    elif ifile2 in ['', None, 'None']:
        merged_df = auto_decode(ifile1)
    else:
        merged_df = pd.merge(auto_decode(ifile1),auto_decode(ifile2),how='outer',left_index=True,right_index=True).fillna(0)
    merged_df = delete_zero_row(merged_df)
    merged_df.index.name = 'Taxonomy'
    return merged_df





# convert mpa to newick tree
def tree(): return defaultdict(tree)
def dicts(tree):
    return {key: (dicts(tree[key]) if hasattr(tree[key], 'items') else tree[key]) for key in tree}

def newickify(node_to_children, root_node, dist) -> str:
    visited_nodes = set()

    def newick_render_node(name, distance: float) -> str:
        assert name not in visited_nodes, "Error: The tree may not be circular!"

        if name not in node_to_children:
            # Leafs
            return F'{name}:{distance}'
        else:
            # Nodes
            visited_nodes.add(name)
            children = node_to_children[name]
            children_strings = [newick_render_node(child, children[child]) for child in children.keys()]
            children_strings = ",".join(children_strings)
            return F'({children_strings}){name}:{distance}'

    newick_string = newick_render_node(root_node, dist) + ';'

    # Ensure no entries in the dictionary are left unused.
    assert visited_nodes == set(node_to_children.keys()), "Error: some nodes aren't in the tree"

    return newick_string
def map2newick(dat, dist):
    tree_hash = tree()
    for i in dat.index:
        tmp = i.split('|')
        tree_hash['root'][tmp[0]] = dist
        for n in range(1,len(tmp)):
            tree_hash[tmp[n-1]][tmp[n]] = dist
        
    newick_tree = newickify(tree_hash, 'root', dist)
    return newick_tree

def filter_metadata_att(metadata,group_id='phenotype'):
    data = auto_decode(metadata)
    if group_id in data.columns:
        if len(data[group_id].unique()) != 2:
            print('Error: the column {} does not have exact 2 level.'.format(group_id))
            exit(1)
    else:
        print('Error: the column {} does not exists in the metadata input.'.format(group_id))
        exit(1)
    
    outputlist=[]
    for att in data.columns:
        if data[att].dtypes == object:
            l = len(data[att].unique())
            if l != 1 and l <= 10:
                outputlist.append(att)
            else:
                print('Warning: Discard metadata attribute {} for the level number of this attribute equal to 1 or greater than 10.'.format(att))
        else:
            if len(data[att].unique()) != 1:
                outputlist.append(att)
        if len(outputlist) == 10:
            print('Warning: Too many attributes in metadata, reserved the first 10 attributes in the output.')
            break
    tmp = data[outputlist]
    return tmp

