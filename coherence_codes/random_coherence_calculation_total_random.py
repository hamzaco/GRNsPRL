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
       l.append((len(gene.is_influenced),len(gene.is_influenced_externals)))
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
with open('coherences_new.json', 'r') as fp:
    coherences=json.load(fp)
entries = os.listdir(path)
error_bars={}
p_values={}
coherences_random={}
for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "/" + "expressions.ALL.txt"
    externals = path + "/" + cell + "/" + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    inputs=NumberOfInputs(original.genes)
    original = Logical_Table_Network(cell, expressions, externals)
    original.readData()
    random_size=100
    coherences_data=[]
    p_values[cell]=0
    for i in range(0,random_size):
        gn=random_network(original,method="brute force")
        gn.calculate_attractors()
        gn.coherence_network()
        #print_genes(gn.genes)
        coherences_data.append(gn.coherence)
        if coherences[cell]<gn.coherence:
            p_values[cell] += 1
        print(gn.coherence)
    coherences_random[cell]=statistics.mean(coherences_data)
    error_bars[cell]=statistics.stdev(coherences_data)

with open('p_values_totally_random_coherences.json', 'w') as fp:
    json.dump(p_values, fp)
with open('coherences_totally_random.json', 'w') as fp:
    json.dump(coherences_random, fp)
with open('error_bars_totally_random_coherences.json', 'w') as fp:
    json.dump(error_bars, fp)
