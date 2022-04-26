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
entries = os.listdir(path_cells)


sens=[]
bias=[]
deviations_random=[]
sample_size=100
count=0
ncf_count=0
count_one_input=0
count_large=0
def boundary_edges(n,h):

    if n==0:
        return 0
    if n==1:
        return 1
    if h%2==0:
        return 2*boundary_edges(n-1,h/2)
    if h>2**(n-1):
        return boundary_edges(n-1,2**n-h)+2**n-h
    else:
        return boundary_edges(n-1,h)+h
deviation_from_one=[]
deviations=[]
"""for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for gene in original.genes:
        count += 1
        if gene.number_of_inputs() >13:
            print(gene.number_of_inputs(), gene)
            count_large += 1
            continue

        if gene.number_of_inputs()<3:
            count_one_input +=1
            continue
        s,input_number=gene.find_sensitivity()
        n=gene.number_of_inputs()
        h = int(gene.activity_bias() * 2 ** n)
        ncf_s=boundary_edges(n,h)/(2**(n-1))
        if ncf_s==0:
            continue
        if ncf_s !=s:
            ncf_count += 1
            deviation_from_one.append(abs(s-1))
        deviations.append((s-ncf_s)/ncf_s)
        for j in range(sample_size):
            gene.logical_dict = random_function(n, number_of_truths=h)
            s = gene.find_sensitivity()[0]
            sens.append(s)
            bias.append(gene.activity_bias())
            deviations_random.append((s - ncf_s) / ncf_s)

print(statistics.mean(deviation_from_one))
print(count,"ncf:",ncf_count)
print("1-input", count_one_input)
print("Large ",count_large)
"""

with open("deviations_biological.txt", "r") as fp:
    deviations=json.load(fp)

with open("deviations_random.txt", "r") as fp:
    deviations_random=json.load(fp)
plt.figure(figsize=(6.4, 4.8))
ax=plt.axes()
bins=ax.hist(deviations,bins=25,range=[0,2.5],log=True
             ,
             histtype="step",alpha=1,label="Cell Collective",linewidth=2)
weights=[1/100]*(len(deviations_random))
bins_random=ax.hist(deviations_random,bins=25, range=[0, 2.5],weights=weights,
         histtype="step",ec="k",log=True,alpha=1,label="Randomized",linewidth=2)
#print(bins[2][0])
#print(bins_random[2][0])
ax.set_xlabel("Normalized Excess Sensitivity "+r'($\delta$)',fontsize=20)
ax.set_ylabel("Count",fontsize=20)
ax.set_xticks(np.arange(0,2.51,step=0.5))
ax.set_xticks(np.arange(0,2.5,step=0.1),minor="true")
ax.tick_params(labelsize=17)
ax.legend(fontsize=22)
ax.set_ylim(0.5,2500)
ax.set_title("\n")
plt.text(1,5000,"(a)",fontsize=20)
plt.tight_layout()
plt.savefig("normalized_excess_ncf_absolute_log_random.png")
plt.show()

ax=plt.axes()
ax.hist(deviations,bins=20,range=[0,1],ec="k",alpha=0.2,label="Cell Collective",color="red")
ax.hist(deviations_random,bins=50, range=[0, 2.5],weights=[1/100]*len(deviations_random),
         ec="k",color="blue",alpha=0.2,label="Random")

ax.set_xlabel("Norm. Excess Sensitivity"+r'($\delta$)',fontsize=20)
ax.set_ylabel("Count",fontsize=20)
ax.set_xticks(np.arange(0,2.51,step=0.5))
ax.set_xticks(np.arange(0,2.5,step=0.1),minor="true")
ax.tick_params(labelsize=15)
ax.set_ylim(1,1000)
ax.legend(fontsize=20)

plt.tight_layout()
plt.savefig("normalized_excess_ncf_absolute_random.png", bbox_inches="tight")
#plt.show()

non_zero=[]
for x in deviations:
    if x!= 0:
        non_zero.append(x)
print(statistics.mean(non_zero))

non_zero=[]
for x in deviations_random:
    if x!= 0:
        non_zero.append(x)
print(statistics.mean(non_zero))
with open("deviations_biological_2excluded.txt", "w") as fp:
    json.dump(deviations, fp)

with open("deviations_random_2excluded.txt", "w") as fp:
    json.dump(deviations_random, fp)