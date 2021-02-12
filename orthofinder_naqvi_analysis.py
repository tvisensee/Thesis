"""
Created on Tue Feb  9 15:11:29 2021
OrthoFinder first analysis
@author: tvise
"""
import sys
import csv

def main():
    
    H_orthogroups, HOG_species_dict = read_HOGs("N0_HOGs_1_31.tsv")
    
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
                current_group = set([]) #empty set to be filled with ensembl ids
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
                            current_group.add(ID)
                HOG_species_dict[HOG_num] = species_set
                H_orthogroups[HOG_num] = [current_group]   
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

main()