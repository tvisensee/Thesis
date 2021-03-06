# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 18:59:26 2020

@author: tvise
"""
import sys

def main():
    
    read_OGs('all_but_dog_orthovenn_clusters.txt')
    
    return

def read_OGs(filename):
    """
    Reads a file where each line is an orthogroup made up of ensembl ids
    "MacacaFascicularis|ENSMFAP00000021776.1	HomoSapiens|ENSP00000456585.1", e.g.
    
    """
    #Read the file
    import csv
    orthogroups = {} #list of sets of orthogroup protein IDs
    OG_num = 0
    OG_species_dict = {}
    
    #This chunk of code opens the tab-delimited file with csv reader & reads it
    with open(filename, newline = '') as f:
        content = csv.reader(f, delimiter = '\t')
        for line in content:
            species_set = set([])
            OG_num = OG_num + 1
            current_group = set([]) #empty set to be filled with ensembl ids
            for item in line:   #for each ensembl ID
                if len(item.strip().split('|')) == 2:
                    species, ID = item.strip().split('|')   #split into species name and ID
                else:
                    continue
                current_group.add(ID)
                species_set.add(species)
            OG_species_dict[OG_num] = species_set
            orthogroups[OG_num] = [current_group]   
    f.close()
    print(OG_num)
    print('Read', filename, 'of', len(orthogroups), 'orthogroups.')
    print('Created list of orthogroups and dictionary of species')

    four_count = 0
    for  value in OG_species_dict.values():
        if len(value) == 4:
            four_count = four_count + 1
    print('there are', four_count, 'groups with all four species.')
    
    return orthogroups, OG_species_dict

main()