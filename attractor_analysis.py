from data_structures.Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures import Gene
import itertools
from matplotlib import pyplot as plt
import os
from data_structures.BooleanNetwork import *

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",
"CC_ErbB_1_4_Receptor_Signaling",
"CC_Immune_System_Model"
]
entries = os.listdir(path_cells)

"""all=0
entropy=[]
sens=[]
ave_in_degree=[]
f=open("entropies.txt","r")
for line in f.readlines():
    entropy.append(float(line.split("\t")[1]))
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    bn=gene_network_to_boolean(original)
    #bn.find_attractors(1000)
    in_degrees=[gene.in_degree for gene in original.genes]
    ave_in_degree.append(sum(in_degrees)/len(in_degrees))
    #entropy.append(bn.get_basin_entropy()),
    #sens.append(bn.get_network_sensitivity())
    #f.write(cell+"\t")
    #f.write(str(entropy[-1])+"\n")
f.close()
plt.plot(ave_in_degree,entropy,"*")
plt.show()"""
entropy=[]
for i in range(4,25):
    t=[]
    print(i)
    for j in range(10):
        bn = BooleanNetwork(i,[2]*i,"ncf")
        bn.find_attractors(1000)
        t.append(bn.get_basin_entropy())
    entropy.append(sum(t)/len(t))

plt.plot(np.arange(4,25),entropy,"*")
plt.show()