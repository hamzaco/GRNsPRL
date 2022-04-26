
from data_structures.Logic_Table_Network import compare_logic_dict, restrict
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
non_monotone=0
non_monotones=[]
non_symmetric=[]
xors=[]
conflictings=[]
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        logic_array=gene.construct_logical_array()
        input_size=int(np.log2(logic_array.size))
        monotonicity=[]
        all +=1
        if input_size<3:
            continue
        xor_count,conflicting=gene.behaviours()
        is_symmetric = [1]*input_size == gene.find_symmetries()[0]
        if is_symmetric:
            non_symmetric.append((gene.name,cell))
        if xor_count>0:
            xors.append((gene.name,cell))
        if conflicting:
            conflictings.append((gene.name,cell))
        for i in range(0,input_size):

            f1=restrict(logic_array,i,True)
            f0=restrict(logic_array,i,False)
            c=compare_logic_dict(f1,f0)
            if c==[False,False]:
                non_monotone += 1
                non_monotones.append((gene.name,cell))
                break
            elif c==[True,False]:
                monotonicity.append("-")
            else:
                monotonicity.append("+")

print(all)
print("Nonmonotonicity percentage: ",non_monotone/all,"Nonmonotone function count: ",non_monotone)

print("No symmetry: ", non_symmetric)
print("No monotonone: ", non_monotones)
print("Xor like: ", xors)
print("Conflicting: ", conflictings)


