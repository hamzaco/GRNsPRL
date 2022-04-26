from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from Gene_Network import Gene_Network, shuffle_network,print_genes
import itertools
def NumberOfInputs(Genes):
    l = []
    for gene in Genes:
       l.append(len(gene.is_influenced)+len(gene.is_influenced_externals))
    return l
path = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
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

entries = os.listdir(path)
inputs_and_sensitivities={}
with open('inputs_and_sensitivities.json','r') as fp:
    inputs_and_sensitivities=json.load(fp)

for i in range(0, len(entries)):

    cell = entries[i]
    #print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "/" + "expressions.ALL.txt"
    externals = path + "/" + cell + "/" + "external_components.ALL.txt"
    grn = Gene_Network(cell, expressions, externals)
    #grn.initialize_effects()
    #s,node_sensitivity=grn.sensitivity_test()
    #l=NumberOfInputs(grn.genes)
    #inputs_and_sensitivities[cell]=(node_sensitivity,l)

plt.figure()
for cell in inputs_and_sensitivities:
    temp=inputs_and_sensitivities[cell]
    ns=temp[0]
    inputs=temp[1]
    for i in range(0,len(inputs)):
        if inputs[i]==3:
            if ns[i]==1.25:
                print(cell)
                print(i+1)

    plt.plot(inputs,ns,'*')

plt.show()
