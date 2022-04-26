from data_structures.Logic_Table_Network import logical_gene,Logical_Table_Network
import json
import os
import numpy as np
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from data_structures.Gene import is_symmetric
import itertools


path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]

entries = os.listdir(path_cells)
xor=0
conflicting=0
all=0


for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        if len(gene.is_influenced)+len(gene.is_influenced_externals)>2:
            xor_count,conflicting_count=gene.behaviours()
            if xor_count >0:
                xor +=1
            if conflicting_count > 0:
                conflicting += 1
            all +=1

print(all)
print("Xor behaviour: ",xor/all)
print("Confliciting behaviour: ",conflicting/all)


