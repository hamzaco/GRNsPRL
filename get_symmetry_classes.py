
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
f=open("symmetry_classes_data.txt","w")
not_ncf=[]
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        _,classes=gene.find_symmetries()
        classes=[c[0] for c in classes]
        for i in range(len(classes)):
            line = cell + "\t" + gene.name
            for g in classes[i]:
                line = line+"\t"+g
            line = line +"\n"
            f.write(line)
f.close()