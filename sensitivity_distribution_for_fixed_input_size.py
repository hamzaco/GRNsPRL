
from data_structures.Gene_Network import Gene_Network
from data_structures.Logic_Table_Network import Logical_Table_Network, logical_gene, random_network
from data_structures.Gene import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import math
def NumberOfInputs(Genes):
    l = []
    for gene in Genes:
       l.append((len(gene.is_influenced),len(gene.is_influenced_externals)))
    return l
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
path = "C:/Users/Lenovo/PycharmProjects/INDEP"
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
"CC_Yeast_Apoptosis"]
entries = os.listdir(path_cells)
with open('sensitivity_data/sensitivity_distribution_original.json','r') as fp:
    original_sensitivities=json.load(fp)

input_numbers=[3]
d=[0,0]


count=[0]*20
count_random=[0]*20
all=[0]*20
"""for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for i in range(0,len(original.genes)):
        gene= original.genes[i]
        s,real_size=gene.find_sensitivity()
        temp=original_sensitivities[real_size].copy()
        temp.append(s)
        original_sensitivities[real_size]=temp
    random_size=1
    gn = Logical_Table_Network(cell, expressions, externals)
    gn.readData()
    inputs=NumberOfInputs(original.genes)
    for j in range(0, random_size):
        gn = random_network(gn, method="symmetry dictated")
        for i in range(0,len(gn.genes)):
            gene= gn.genes[i]
            s,real_size = gene.find_sensitivity()

            temp = random_sensitivities[real_size].copy()
            temp.append(s)
            random_sensitivities[real_size] = temp"""

for input_size in range(1,len(original_sensitivities)):
    plt.figure()
    originals=original_sensitivities[input_size]
    if len(originals)== 0:
        continue
    plt.hist(originals, bins=40, range=[0,4],weights=[1/len(originals)]*len(originals), facecolor='red', ec='black', alpha=0.2, label='Original')
    plt.legend()
    plt.title(str(input_size)+' - input nodes sensitivity distribution')
    plt.xlabel('Sensitivity coefficient')
    plt.ylabel('Probability')
    plt.show()


