from data_structures.Logic_Table_Network import Logical_Table_Network
import json
import os
import numpy as np
import scipy.optimize as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, Gene
from scipy.stats import pearsonr
import random
import itertools
from sample_functions import random_function

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",]
entries = os.listdir(path_cells)

sens_activity_bias=[]
activity_bias=[]
plt.figure()
canalizing_layers=[[]]*20

def boundary_edges(n,h):
    if n==1:
        return 1
    if h>2**(n-1):
        return boundary_edges(n-1,2**n-h)+2**n-h
    else:
        return boundary_edges(n-1,h)+h

for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for gene in original.genes:
        if gene.number_of_inputs()==1:
            continue
        canalizing_nodes,_=gene.is_canalizing()
        for j in range(len(canalizing_nodes)-1):
            temp=canalizing_layers[j+1].copy()
            temp.append(canalizing_nodes[j][2])
            canalizing_layers[j+1]=temp

observed_prob=[0]*9
for i in range(9):
    if len(canalizing_layers[i]) != 0:
        observed_prob[i]=sum(canalizing_layers[i])/len(canalizing_layers[i])

#alpha=sp.minimize_scalar(lambda a: sum(abs(np.array(observed_prob)-
 #                             np.array([2*np.exp(-2**(-m)*a)/(1+np.exp(-2**(-m)*a)) for m in range(14)])))).x

alpha=sp.minimize(lambda a: sum(abs(np.array(observed_prob)-
                              np.array([a[1]*np.exp(-2**(-m)*a[0])/(1+np.exp(-2**(-m)*a[0])) for m in range(9)]))),[1,1]).x

guess_prob= [alpha[1]*np.exp(-2**(-m)*alpha[0])/(1+np.exp(-2**(-m)*alpha[0])) for m in range(9)]
print(alpha)
yerror=[1/np.sqrt(4*len(x)) for x in canalizing_layers[1:9]]
plt.plot(np.arange(1,9),observed_prob[1:],"*k",label="Observed")
plt.errorbar(np.arange(1,9),np.array(observed_prob[1:]),yerr=np.array(yerror),fmt="*k")
plt.plot(np.arange(1,9),guess_prob[1:],"*r",label="Guessed Expression")
plt.xticks(np.arange(0,9,step=1))
plt.legend()
plt.ylabel("Probability")
plt.xlabel("Layer ")
plt.savefig("ncf_layer_probability_fit.png")
plt.show()