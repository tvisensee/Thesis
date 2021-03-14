# -*- coding: utf-8 -*-
"""
Created on Sun Mar  7 19:53:26 2021

@author: tvise
"""
#import orthofinder_naqvi_analysis as orfind
import csv

def main():
    
    name_dict = read_ID_key("sexbias.genename.ensembl.IDs.csv")
    
    ID_dict, gene_dict, info = read_mappings("HOG_mappings.csv")
    
    map_dict = find_unmapped_IDs(name_dict, ID_dict, gene_dict)
    
    return

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
                #print(header)
                first_line = False
            elif first_line == False:
                for i in range(len(line)):
                    name = line[0]
                    ensembl_IDs = line[1:len(line)]
                    name_dict[name] = ensembl_IDs
    f.close()
    ID_gene_dict = {}
    for gene, OG in name_dict.items():
        for ID in OG:
            ID_gene_dict[ID] = gene
            
    print("\nread_ID_key() on", filename, "making translation dictionary with", len(name_dict), "genes & IDs")
    return ID_gene_dict

def read_mappings(filename):
    ID_dict = {}
    gene_dict = {}
    info = {}
    print('\nread_mappings() on', filename, 'commencing')
    first_line = 'headers'
    with open(filename, newline = '') as f:
        content = csv.reader(f, delimiter = ',')
        for line in content:
            if first_line == 'headers':
                first_line = 'done'
            elif first_line == 'done':
                gene = line[0]
                OG_num = line[1]
                species_num = line[2]
                tissue = line[3]
                bias = line[4]
                OG_contents = line[6:len(line)]
            
                ID_dict[OG_num] = OG_contents
                gene_dict[gene] = OG_contents
                info[OG_num] = [gene, species_num, tissue, bias]
    print('\nRead in', len(gene_dict), 'orthogroups from', filename)            
    return ID_dict, gene_dict, info

# function finds ids within orthogroups that do not map to the corresponding
# Naqvi gene; outputs a dictionary where keys are genes and values are lists
# of [mapped_IDs, unmapped_IDs]
def find_unmapped_IDs(name_dict, ID_dict, gene_dict):
    print('\nfind_unmapped_IDs() commencing')
    map_dict = {}
    unmapped_count = 0
    for gene, OG in gene_dict.items():
        mapped_IDs = []
        unmapped_IDs = []
        for ID in OG:
            if ID in name_dict.keys():
                mapped_gene = name_dict[ID]
                if gene == mapped_gene:
                    mapped_IDs = mapped_IDs + [ID]
                else:
                    unmapped_IDs = unmapped_IDs + [ID]
            else:
                unmapped_IDs = unmapped_IDs + [ID]
                unmapped_count = unmapped_count + 1
        map_dict[gene] = [mapped_IDs, unmapped_IDs]
    print(map_dict)
    print('\nFound', unmapped_count, 'unmpapped IDs.')
    print('\nCreated a', len(map_dict.keys()),'entry dictionary of mapped and unmapped IDs')
    return map_dict
'''
def merge_dicts(map_dict, info_dict):
    # final output structure = dict[OG#] = [gene, [ids that map to gene], [IDs that don't map to gene]]
    merged = {}
    ID_OG = {}
    
    for OG_num, contents in info_dict.items():
        for ID in contents:
            ID_OG[ID] = OG_num
    
    for 
    
    return

def write_mappings_readable(filename, dictionary):
    
    return

def write_mappings_tidy():
    
    return

'''
main()