
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

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= ["CC_CD4p_Tcell_Diff",
"CC_Differentiation_T_lymphocytes",
"CC_Inflammatory_Bowel_Disease_IBD_Model",
"CC_Glucose_Repression_Signaling_2009",
"CC_HCC1954_Breast_Cell_Line_Long-term_ErbB_Network",
"CC_IGVH_mutations",
"CC_Lymphopoiesis_Regulatory_Network",
"CC_Signaling_in_Macrophage_Activation",
"CC_Signal_Transduction_in_Fibroblasts",
"CC_Treatment_of_Castration_Resistant_Prostate_Cancer",
              "CC_EGFR_ErbB_Signaling",
"CC_Yeast_Apoptosis"]

entries = os.listdir(path_cells)
canalizing_percentages=[]
s_vs_canalizing=[]
all=0
c=0
for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    canalizing_percentage=0
    Genes = original.genes
    canalizing_count = 0
    gn = Logical_Table_Network(cell, expressions, externals)
    gn.readData()
    Genes=gn.genes
    for i in range(len(Genes)):
        gene=Genes[i]
        o_gene= original.genes[i]
        gene.is_influenced= o_gene.is_influenced
        gene.is_influenced_externals = o_gene.is_influenced_externals
        #gene.logical_dict=gene.construct_logical_dict()
        #dist,_=gene.find_symmetries()
        #s = gene.find_sensitivity()
        if not gene.is_canalizing():
            c += 1
            canalizing_count +=1
        all += 1
    #canalizing_percentage=canalizing_count/len(Genes)
    #canalizing_percentages.append(canalizing_percentage)
    #s_vs_canalizing.append((canalizing_percentage,s))
#plt.hist(canalizing_percentages,bins=10,range=[0,1],facecolor='orange', ec='black',weights=[1/len(canalizing_percentages)]*len(canalizing_percentages))
#plt.xlabel('Canalizing percentage in a network')
#plt.ylabel('Normalized frequency')
#plt.savefig('canalizing_percentages_distribution')
#plt.show()
print(c)