"""
Created on Wed Oct 21 16:03:25 2020
##
Get tissue-specific sex-biased genes from Naqvi et al's (2019) sex_bias_conserved.matrix.txt
and write to a new file 
##
@author: tvise
"""
def main():
    
    tis_gene_exp = read_matrix('sexbias.conserved.matrix.txt')
    
    tissues = ['Skin', 'Heart', 'Brain', 'Muscle', 'Spleen', 'Adipose', 'Thyroid', 'Colon', 'Lung', 'Liver', 'Pituitary', 'Adrenal']
    
    #write_info(tis_gene_exp, '.sexbias.conserved.genes.txt', tissues, genesOnly=True)
    stats_table = get_stats(tis_gene_exp)
    
    write_file(stats_table, 'sbc_quicksum.csv')
    
    return

#copied and adapted read_matrix() from part2.py from bio331 
def read_matrix(filename):
    """
    Reads a file provided as a tab-delimited matrix of tissues (columns), genes (rows),
    and expression biases as a range from 
    """
    #Read the file and populate nodes (list) and an adjacency matrix 
    #(list of lists) with values representing contact time between badgers
    import csv
    gene_mat = []
    
    #This chunk of code opens the tab-delimited file with csv reader & reads it
    with open(filename, newline = '') as f:
        matrix = csv.reader(f, delimiter = '\t')
        first_line = True
        gene_dict = {} #tissue: [genes expressed in tissue: expression value]
        for line in matrix:
            if first_line == True:
                tissue = line[0:len(line)]
                print(tissue)
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
                            gene_dict[tis_dict[i]] = gene_dict[tis_dict[i]] + [exp_list]
    for key in gene_dict.keys():
        print('In',key, 'there are', len(gene_dict[key]), 'conserved sex-biased genes.')
    f.close()
    print('\nSex-biased genes are now in a dictionary of {tissue: [[gene, expression bias]]}')
    return gene_dict

def get_stats(dictionary):
    #prints stats about matrix
    stats_table = [['Tissue', 'GeneCount', 'FemaleBias', 'MaleBias']]
    print('\n' + 'Tissue' + '     ' + 'Genes' + '     ' +'Fb_genes' + '     ' + 'Mb_genes', end = '')
    for tissue in dictionary.keys():
        female = 0
        male = 0
        for gene in dictionary[tissue]:
            count = len(dictionary[tissue])
            if gene[1] < 0:
                female = female + 1
            elif gene[1] > 0:
                male = male + 1
        row = [tissue,count,female,male]
        stats_table = stats_table + [row]
        print("\n" + tissue + "     " + str(count) + '     ' + str(female) + '     ' + str(male), end = "")
    
    return stats_table

def write_file(table, outfile):
    import csv
    with open(outfile, 'w', newline = '') as f:
        info = csv.writer(f, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)
        info.writerows(table)
        
    print('\nwrote table to', outfile)
    return

def write_info(dictionary, outfile, tissue, genesOnly=None):
    """
    writes a tab-delimited .txt file from gene_dict of gene names and expression values
    outputs different files for specific tissues
    """
    import csv
    #for every tissue in list of supplied tissue, generate an output file of genes and bias scores
    for item in tissue:
        #Write \t delimited file with csv writer
        name = item + outfile
        with open(name, 'w', newline = '') as f:
            info = csv.writer(f, delimiter = '\t', quoting=csv.QUOTE_NONE, skipinitialspace = True) #creates object to be written to
            if genesOnly == True:
                genes = []
                for gene in dictionary[item]:
                    genes = genes + [gene[0]]
                info.writerow(genes) # write all gene names on seperate rows
                print("\nWrote only genes to", name) 
            elif genesOnly == False:
                info.writerows(dictionary[item]) #write the [gene, bias] lists as rows to the file
                print("\nWrote genes and bias scores to", name) 
    return 

if __name__ == '__main__':
    main()