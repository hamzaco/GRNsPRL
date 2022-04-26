from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene import Gene
from Gene_Network import Gene_Network, shuffle_network,print_genes
from scipy.stats import pearsonr
import itertools
x=Gene('A')
y=Gene('B')
z=Gene('C')
t=Gene('D')
w=Gene('E')

g1= Gene(name='test1',function='A and B and not C' )
g1.is_influenced=set([x,y,z])

Genes=[x,y,z,t,w,g1]
dist,_=g1.find_symmetries()
print(dist)

"""
def NumberOfInputs(Genes):
    l = []
    for gene in Genes:
       l.append((len(gene.is_influenced),len(gene.is_influenced_externals)))
    return l
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
sensitivities_random_symmetric_wo_xor={}
partitions_random_wo_xor=[[]]*20
activity_bias_random=[]
with open('sensitivities.json', 'r') as fp:
    sensitivities = json.load(fp)

with open('activity_bias.json','r') as fp:
    activity_bias_original=json.load(fp)

dist_vs_input_number_original=[[]]*20
conf_number=0
all=0
summ=0
all_count=0
count_sens=0
count_str=0

plt.figure()
for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    Genes = original.genes
    inputs = NumberOfInputs(original.genes)

    for gene in Genes:
        s, realsize = gene.find_sensitivity()
        dist,classes=gene.find_symmetries()
        func, op_rel = gene.find_structure()
        print(classes)
        wait=input('say:  ')
        bias = gene.number_of_truths()
        if realsize <= 2:
            continue
        s_p=0
        for i in range(len(dist)):
            temp= 2**(realsize-dist[i])*2*bias*(1-bias)
            s_p += temp*dist[i]/(2**(realsize-1))
        plt.plot(bias,s,'*k')
        plt.plot(bias,s_p,'*r')
        all_count +=1
        bias=gene.number_of_truths()
        s_m=gene.my_sensitivity()

        if s!=s_m:
            count_sens+=1
            print(s_m)
            print(s)

        func, op_rel= gene.find_structure()

        if func=='not found':
            count_str += 1

plt.show()
with open('dist_vs_input_number_original.json', 'w') as fp:
    json.dump(dist_vs_input_number_original, fp)
"""
