"""
Created on Tue Feb  9 15:11:29 2021
OrthoFinder first analysis
@author: tvise
"""
import sys
import csv

def main():
    
    H_orthogroups, HOG_species_dict = read_HOGs("N0_HOGs_1_31.tsv")
    
    gene_dict, gene_list = read_matrix("sexbias.conserved.matrix.txt")
    
    name_dict = read_ID_key("sexbias.genename.ensembl.IDs.csv")
    
    map_gene_names(gene_list, name_dict, H_orthogroups)
    return

#copied from orthovenn_analysis.py and adapted to orthofinder files
def read_HOGs(filename):
    """
    Reads a file where each line is an orthogroup
    Column headers are:
        
    HOG	OG	Gene Tree Parent Clade	Canislupusfamiliaris	Homosapiens	
    Macacafascicularis	Musmusculus	Rattusnorvegicus
    
    If there are more than 1 genes from a single species they will be comma-
    separated under the species column name 
    """
    #Read the file
    H_orthogroups = {}
    H_orthogroups_rev = {}
    HOG_species_dict = {}
    first_line = "headers"
    HOG_num = 0
    
    #This chunk of code opens the tab-delimited file with csv reader & reads it
    with open(filename, newline = '') as f:
        content = csv.reader(f, delimiter = '\t')
        for line in content:
            if first_line == "headers":
                first_line = "done"
                continue
            elif first_line == "done":
                species_set = set([])
                current_group = [] #empty list to be filled with ensembl ids
                #define important columns
                HOG_num = line[0]
                index = 4 #index/column at which the species begin
                species_list = ["dog", "human", "macaque", "mouse", "rat"]
                item_count = 0
                for item in line:   #for each entry
                    item_count = item_count + 1
                    if item_count >= index and len(item) > 0:
                        species_set.add(species_list[item_count - index])
                        IDs = item.strip().split(',')
                        for ID in IDs:
                            trunc_ID = truncate_ID(ID) #strips off "." and anything following
                            current_group = current_group + [trunc_ID]
                HOG_species_dict[HOG_num] = species_set
                H_orthogroups[HOG_num] = current_group
    f.close()
    
    print(HOG_num)
    print('Read', filename, 'of', len(H_orthogroups), 'orthogroups.')
    print('Created list of orthogroups and dictionary of species')

    five_count = 0
    for  value in HOG_species_dict.values():
        if len(value) == 5:
            five_count = five_count + 1 
    print('there are', five_count, 'groups with all five species.')
    print('Example hierarchical orthogroup dictionary entry:\n', "N0.HOG0015353:\n",
          H_orthogroups["N0.HOG0015353"])
    
    return H_orthogroups, HOG_species_dict

#copied and adapted read_matrix() from orthovenn_analysis.py 
def read_matrix(filename):
    """
    Reads a file provided as a tab-delimited matrix of tissues (columns), genes (rows),
    and expression biases as a range from -1 to 0 (female biases are negative)
    """
    #This chunk of code opens the comma-delimited file with csv reader & reads it
    with open(filename, newline = '') as f:
        matrix = csv.reader(f, delimiter = '\t')
        first_line = True
        gene_dict = {} #tissue: [genes expressed in tissue: expression value]
        gene_list = [] #list of sex-biased conserved genes
        for line in matrix:
            if first_line == True:
                tissue = line[0:len(line)]
                tis_dict = {} # index(int): tissue(str) 
                for i in range(len(tissue)):
                    tis_dict[i + 1] = tissue[i]
                gene_dict = {t: [] for t in tissue}
                first_line = False
            #then each line in the file is
            elif first_line == False:
                exp_list = [] # expression list, gene(str), bias (int)
                for i in range(len(line)):
                    gene = line[0]
                    if i > 0:
                        bias = int(line[i])
                        if bias == 0:
                            pass
                        elif bias != 0:
                            exp_list = [gene, bias]
                            gene_list = gene_list + [gene]
                            gene_dict[tis_dict[i]] = gene_dict[tis_dict[i]] + [exp_list]
    f.close()
    for key in gene_dict.keys():
        print('In',key, 'there are', len(gene_dict[key]), 'conserved sex-biased genes.')
    print("\nIn total, there are", len(gene_list), "sex-biased, conserved genes.")
    return gene_dict, gene_list

def read_ID_key(filename):
    # Read file of common gene names and the 5 species ensembl gene IDs into 
    # a dictionary
    with open(filename, newline = '') as f:
        content = csv.reader(f, delimiter = ',')
        first_line = True
        name_dict = {} # dictionary of gene name: [list of ensembl IDs]
        for line in content:
            if first_line == True:
                header = line[0:len(line)]
                print(header)
                first_line = False
            elif first_line == False:
                for i in range(len(line)):
                    name = line[0]
                    ensembl_IDs = line[1:len(line)]
                    name_dict[name] = ensembl_IDs
    f.close()
    print("Read file with", len(name_dict), "common gene names and associated ensembl protein IDs")
    return name_dict

def map_gene_names(gene_list, name_dict, OGs):
    names_found = 0
    not_found = 0
    mapped_IDs = 0
    OG_hit_set = set([])
    OG_hit_dict = {}
    OG_gene_dict = {}
    for gene in gene_list:
        if gene in name_dict.keys():
            ensembl_IDs = name_dict[gene]
            names_found = names_found + 1 
            for ID in ensembl_IDs:
                if names_found == 1:
                    print(ID)
                else:
                    continue
                for OG_num, group in OGs.items():
                    if ID in group: 
                        mapped_IDs = mapped_IDs + 1
                        OG_hit_set.add(OG_num)
                        OG_gene_dict[OG_num] = [gene, group]
                        print("hit found!", OG_num, gene, len(group))
        else:
            not_found = not_found + 1
            continue
    print(names_found, "sex biased, conserved gene names found in the name translation dict")
    print(not_found, "sbc gene names not found in the name translation dict")
    print("Orthogroups that genes mapped to:", len(OG_hit_set))
    return

def truncate_ID(ID):
    
    trunc_ID, suffix = ID.split(".")
    
    return trunc_ID

#adapted from geeksforgeeks.org
def get_key(val, dictionary):
    for key, value in dictionary.items():
        if val == value:
            return key
        else:
            return "no key found for value"

main()