
from data_structures.Logic_Table_Network import compare_logic_dict, restrict, random_network
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

not_ncf=[]
asymmetric=[]
not_canalizing=[]
all=0
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        if gene.number_of_inputs()<3:
            continue
        all += 1
        nodes, is_ncf=gene.is_canalizing()
        if not gene.is_symmetric():
            asymmetric.append(gene)
        if not is_ncf:
            not_ncf.append((gene.name,cell)),
        if len(nodes)==0:
            not_canalizing.append((gene.name,cell))
        else:
            s, _ = gene.find_sensitivity()
            a=gene.activity_bias()
            plt.plot(a,s,'*k')

plt.show()
print(len(not_ncf)," not nested canalizing out of",all)
print(len(asymmetric),"asymmetric out of", all)
print(len(not_canalizing),"not canalizing out of",all)
