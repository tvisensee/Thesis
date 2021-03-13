"""
Created on Tue Feb  9 15:11:29 2021
OrthoFinder first analysis
@author: tvise
"""
import sys
import csv

def main():
    
    H_orthogroups, HOG_species_dict = read_HOGs("N0_HOGs_1_31.tsv")
    
    gene_tissue, gene_list, exp_dict, no_tissue_genes = read_matrix("sexbias.conserved.matrix.txt")
    
    name_dict = read_ID_key("sexbias.genename.ensembl.IDs.csv")
    
    OG_hit_set, OG_hit_dict, OG_gene_dict, OG_mapped_IDs = map_gene_names(gene_list, name_dict, H_orthogroups)
    
    file_dict = merge_dicts(OG_hit_dict, OG_gene_dict, gene_tissue, HOG_species_dict, H_orthogroups, exp_dict)
    
    write_dict("HOG_mappings.csv", file_dict)
    
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

    print('\nread_HOGs() on', filename, 'of', len(H_orthogroups), 'orthogroups.')
    print('Created list of', len(H_orthogroups), ' orthogroups and dictionary of species')

    five_count = 0
    for  value in HOG_species_dict.values():
        if len(value) == 5:
            five_count = five_count + 1 
    print('there are', five_count, 'orthogroups with all five species.')
    print('Example hierarchical orthogroup dictionary entry:\n', "N0.HOG0015353:\n",
          H_orthogroups["N0.HOG0015353"])
    
    return H_orthogroups, HOG_species_dict

#copied and adapted read_matrix() from orthovenn_analysis.py 
def read_matrix(filename):
    print("\nread_matrix() on", filename)
    """
    Reads a file provided as a tab-delimited matrix of tissues (columns), genes (rows),
    and expression biases as a single value: 0 (no bias), -1, or 1 (female biases are negative)
    """
    #This chunk of code opens the comma-delimited file with csv reader & reads it
    with open(filename, newline = '') as f:
        matrix = csv.reader(f, delimiter = '\t')
        first_line = True
        gene_tissue = {} #gene: tissue
        exp_dict = {}
        gene_list = [] #list of sex-biased conserved genes
        for line in matrix:
            if first_line == True:
                tissue = line[0:len(line)]
                tis_dict = {} # index(int): tissue(str) 
                for i in range(len(tissue)):
                    tis_dict[i + 1] = tissue[i] #populate tis_dict with tissues/headers
                first_line = False
            #then each line in the file 
            elif first_line == False:
                gene = line[0]
                for i in range(len(line)):
                    if i > 0:
                        bias = int(line[i])
                        if bias == 0:
                            pass
                        elif bias != 0:
                            exp_dict[gene] = bias
                            gene_list = gene_list + [gene]
                            gene_tissue[gene] = tis_dict[i]
    f.close()
    print("Read", filename, "and created list of", len(gene_list), "sex-biased, conserved genes.")
    print("Created a gene: tissue dictionary containing", len(gene_tissue.keys()), "genes with tissues.\n")
    
    no_tissue_genes = []
    for gene in gene_list:
        if gene not in gene_tissue.keys():
            no_tissue_genes = no_tissue_genes = [gene]
    
    return gene_tissue, gene_list, exp_dict, no_tissue_genes

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
    print("\nread_ID_key() on", filename, "making translation dictionary with", len(name_dict), "genes & IDs")
    return name_dict

# function takes information from name dictionaries to convert from common gene
# names to ensembl IDs and then map which OGs contain those ensembl IDs
def map_gene_names(gene_list, name_dict, OGs):
    print('\nmap_gene_names() starting to map', len(gene_list), 'genes to their orthogroups')
    names_found = 0
    not_found = 0
    no_og = 0
    OG_hit_set = set([])
    OG_hit_dict = {}
    OG_gene_dict = {}
    OG_mapped_IDs = {}
    found = False
    #for every conserved sex-biased gene name
    for gene in gene_list:
        # check if gene is in name-to-ensemblID dictionary
        if gene in name_dict.keys():
            ensembl_IDs = name_dict[gene] #translate gene to ensemblIDs (genes may map to multiples species IDs)
            names_found = names_found + 1 
            found = True   # mark gene as found
            # for all ensemblIDs that found gene maps to
            for ID in ensembl_IDs:
                if found == True:
                    #print(ID)
                    # find OG containing those IDs
                    for OG_num, group in OGs.items():
                        mapped_IDs = []
                        unmapped_IDs = group
                        if ID in group: 
                            # store info about where this gene maps
                            OG_hit_set.add(OG_num)
                            OG_hit_dict[gene] = OG_num
                            OG_gene_dict[OG_num] = [gene, group]
                            mapped_IDs = mapped_IDs + [ID]
                            unmapped_IDs = unmapped_IDs.remove(ID)
                        OG_mapped_IDs[OG_num] = [mapped_IDs, unmapped_IDs]
                            #print("hit found!", OG_num, gene, len(group))
                else:
                    # keep track of genes found in name dict that don't map to any OGs
                    no_og = no_og + 1
                    continue
        else:
            # keep track of things not found in name dictionary
            found = False
            not_found = not_found + 1
            continue

    print(names_found, "sex biased, conserved gene names found in the name translation dict")
    print(not_found, "sbc gene names not found in the name translation dict")
    print(no_og, "genes not mapped to any orthogroup")
    print("Orthogroups that genes mapped to:", len(OG_hit_set))
    
    return OG_hit_set, OG_hit_dict, OG_gene_dict, OG_mapped_IDs

#takes info stored in disparate dictionaries and merges them into the format:
    # gene name: [OG#, # species, tissue, OG size, bias, OG content]
# note that naming format for the input dictionaries is key_value(_value)
def merge_dicts(gene_OGnum, OGnum_gene_group, gene_tissue, OGnum_species, OGnum_group, gene_bias):
    merge = {}
    
    for gene in gene_OGnum.keys():
        OGnum = gene_OGnum[gene]
        num_species = len(OGnum_species[OGnum])
        #there may not be a tissue for every gene
        try:
            tissue = gene_tissue[gene]
        except KeyError:
            tissue = 'NA'
        OG_size = len(OGnum_group[OGnum])
        bias = gene_bias[gene]
        OG_content = OGnum_group[OGnum]
        
        info_list = [OGnum, num_species, tissue, OG_size, bias] + OG_content
        merge[gene] = info_list
    
    return merge

def write_dict(filename, dictionary):
    print('\nwrite_dict() writing to', filename)
    filerows = 0
    with open(filename, 'w', newline = '') as f:
        info = csv.writer(f, delimiter = ',', quoting=csv.QUOTE_NONE, skipinitialspace = True)
        row1 = ["Gene", "OG_number", "OG_species", "Tissue", "OG_size", "Bias", "OG_contents"]
        info.writerow(row1)
        for gene in dictionary.keys():
            row = [gene] + dictionary[gene]
            filerows = filerows + 1
            info.writerow(row)
    f.close()
    print('Wrote', filerows, 'rows to', filename)
    return

def write_unmapped(filename, dictionary):
    print('\nwrite_unmapped() writing to', filename)
    filerows = 0
    unmapped_IDs = {}
    for OG in dictionary.keys():
        unmapped = dictionary[OG][1]
        if len(unmapped) > 0:
            unmapped_IDs[OG] = unmapped
    with open(filename, 'w', newline = '') as f:
        info = csv.writer(f, delimiter = ',', quoting=csv.QUOTE_NONE, skipinitialspace = True)
        row1 = ['OG', 'Unmapped_IDs']
        info.writerow(row1)
        for OG in unmapped_IDs.keys():
            row = [OG] + [unmapped_IDs[OG]]
            filerows = filerows + 1
            info.writerow(row)
    f.close()
    print('Wrote', filerows, 'rows to', filename)
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