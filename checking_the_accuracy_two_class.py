
from data_structures.Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from scipy.stats import pearsonr
import itertools

main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"


path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]
types=[0,0,0,0,0]
all=0

plt.figure()
entries = os.listdir(path_cells)

for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        dist, _ = gene.find_symmetries()
        if len(dist) > 2:
            continue
        [s,n]=gene.find_sensitivity()
        if n < 3:
            continue
        all += 1
        a=min(list(dist))
        if s== n/(2**(n-1)):
            types[0] += 1
        elif s==a*2**(1-n)+(n-a)*(2**a-1)*2**(1-n):
            types[1] += 1
        elif s== a*(2**(n-a)-1)*2**(1-n)+(n-a)*2**(1-n):
            types[2] += 1
        elif s== a*(2**(n-a)-1)*2**(1-n)+(n-a)*(2**a-1)*2**(1-n):
            types[3] += 1
        else:
            types[4] += 1
print(all)
print(types)
