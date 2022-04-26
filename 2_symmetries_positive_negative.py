
from data_structures.Logic_Table_Network import compare_logic_dict, restrict, random_network
import json
import os
import numpy as np
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from data_structures.Gene import is_symmetric
import math
import itertools
def nCr(n,r):
    return math.factorial(n)/(math.factorial(n-r)*math.factorial(r))

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]

entries = os.listdir(path_cells)

all2=0
target=0
theoretical=0
upper_estimate_for_variance=0
p_values=[]

for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for gene in original.genes:
        dist,eq_classes=gene.find_symmetries()
        positives=set(gene.positive_signals)
        negatives=set(gene.negative_signals)
        if len(eq_classes)==2:
            all2 += 1
            p=1/nCr(gene.in_degree,dist[0])
            theoretical += p
            p_values.append(p)
            if (positives==set(eq_classes[0]) and negatives==set(eq_classes[1])) or \
                    (positives == set(eq_classes[1]) and negatives == set(eq_classes[0])):
                target += 1


pbar=statistics.mean(p_values)
s=statistics.stdev(p_values)
variance=len(p_values)*pbar*(1-pbar)-len(p_values)*s*s
print(all2)
print(target)
print(theoretical)
print(variance)
n=len(p_values)
T=[0]*n

for i in range(n):
    summ=0
    for j in range(n):
        summ += (p_values[j] / (1 - p_values[j])) ** i
    T[i]=summ

pmf=[0]*n
prod=1
for i in range(n):
    prod *= 1-p_values[i]
pmf[0]=prod
for k in range(1,n):
    summ=0
    for i in range(k):
        summ += (-1)**i*pmf[k-i-1]*T[i+1]
    pmf[k]=summ/k

def pdf(k):
    summ=0
    for i in range(k+1):
        summ += pmf[i]
    return summ

plt.plot(pmf)
plt.savefig("t.png")