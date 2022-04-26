from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from Gene_Network import Gene_Network, shuffle_network,print_genes
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
              "CC_HIV1_interactions_with_TCell_Signalling_Pathway",
"CC_Yeast_Apoptosis"]

entries = os.listdir(path_cells)
sensitivities_canalizing={}
sensitivities_symmetric={}
partitions_canalizing=[[]]*20
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
with open(main_dir+'/sensitivity_data/sensitivities.json', 'r') as fp:
    sensitivities = json.load(fp)
sensitivities.pop("CC_EGFR_ErbB_Signaling")
sensitivities.pop("CC_HIV1_interactions_with_TCell_Signalling_Pathway")

def kullback_leibler(p,q):
    div=0
    for i in range(len(p)):
        div += p[i]*np.log(p[i]/q[i])
    return div

def jensen_shannon(p,q):
    div = 0
    for i in range(len(p)):
        m=(p[i]+q[i])/2
        if p[i]==0 and q[i]==0:
            div += 0
        elif p[i]==0:
            div += (q[i]/2) * np.log(q[i]/m)
        elif q[i]==0:
            div += (p[i]/2) * np.log(p[i] / m)
        else:
            div += (q[i]/2) * np.log(q[i]/m) + (p[i]/2) * np.log(p[i] / m)
    return div

def NumberOfInputs(Genes):
    l = []
    for gene in Genes:
       l.append((len(gene.is_influenced),len(gene.is_influenced_externals)))
    return l


for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    inputs = NumberOfInputs(original.genes)
    gn = Logical_Table_Network(cell, expressions, externals)
    gn.readData()
    random_size = 10
    sensitivities_data = []
    sensitivities_data1 = []

    for i in range(random_size):
        summ = 0
        canalizing= random_network(gn,method="canalizing")
        for gene in canalizing.genes:
            s, real_size = gene.find_sensitivity()
            summ += s
        sensitivities_data.append(summ / len(canalizing.genes))


    for i in range(random_size):
        summ = 0
        symmetric= random_network(gn,method="xor excluded")
        for gene in symmetric.genes:
            s, _ = gene.find_sensitivity()
            summ += s
        sensitivities_data1.append(summ / len(symmetric.genes))
    sensitivities[cell]= sensitivities[cell][0]
    sensitivities_symmetric[cell]=sum(sensitivities_data1)/len(sensitivities_data1)
    sensitivities_canalizing[cell] = sum(sensitivities_data) / len(sensitivities_data)
sensitivities=list(sensitivities.values())
canalizing_sensitivities=list(sensitivities_canalizing.values())
symmetric_sensitivities=list(sensitivities_symmetric.values())

plt.figure()
n_original,_,_=plt.hist(sensitivities,bins=20,range=[0.5, 1.5],weights=[1/len(sensitivities)]*len(sensitivities),
         facecolor='red', ec='black', alpha=0.2, label='Original')
n_canalizing,_,_= plt.hist(canalizing_sensitivities,bins=20,range=[0.5, 1.5],
         weights=[1/len(canalizing_sensitivities)]*len(canalizing_sensitivities),
         facecolor='blue', ec='black', alpha=0.2, label='Canalizing')

plt.legend()
plt.savefig('Sensitivities_histogram_original_canalizing')
plt.show()
plt.figure()
n_original,_,_=plt.hist(sensitivities,bins=20,range=[0.5, 1.5],weights=[1/len(sensitivities)]*len(sensitivities),
         facecolor='red', ec='black', alpha=0.2, label='Original')
n_symmetric,_,_= plt.hist(symmetric_sensitivities,bins=20,range=[0.5, 1.5],
         weights=[1/len(symmetric_sensitivities)]*len(symmetric_sensitivities),
         facecolor='green', ec='black', alpha=0.2, label='Symmetric')
plt.legend()
plt.savefig('Sensitivities_histogram_original_symmetric_xor_excluded')
plt.show()
print("Jensen-Shannon Divergences")
print("Original vs symmetric",jensen_shannon(n_original,n_symmetric))
print("Original vs canalizing",jensen_shannon(n_original,n_canalizing))

ones=[1]*len(canalizing_sensitivities)
error_canalizing=0
error_symmetric=0
error_original=0
for i in range(len(ones)):
    error_canalizing += abs(1-canalizing_sensitivities[i])
    error_symmetric += abs(1-symmetric_sensitivities[i])
    error_original += abs(1-sensitivities[i])

print("Average distances to 1")
print("Symmetric",error_symmetric/len(ones))
print("Canalizing",error_canalizing/len(ones))
print("Original",error_original/len(ones))