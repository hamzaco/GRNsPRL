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
sample_size=10
sens_bias=[]
real_sens_mean=[]
random_sens_mean=[]
random_sens_yerr=[]
"""for i in range(0,len(entries)):
    random_sens_ave=[]
    c=[]
    cell=entries[i]
    print(cell)
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    flag=False
    for j in range(sample_size):
        temp=[]
        for gene in original.genes:
            n=gene.number_of_inputs()
            if n>13:
                continue
            if not flag and not gene.is_canalizing()[1]:
                flag=True
            if j==0:
                c.append(gene.find_sensitivity()[0])
            if n<2:
                temp.append(1)
                continue

            h=int(gene.activity_bias()*(2**n))
            a=Gene()
            a.logical_dict = random_function(n,number_of_truths=h)
            temp.append(a.find_sensitivity()[0])
            sens_bias.append((a.find_sensitivity()[0],h/(2**n)))
        random_sens_ave.append(statistics.mean(temp))
    if flag:
        random_sens_mean.append(statistics.mean(random_sens_ave))
        random_sens_yerr.append(statistics.stdev(random_sens_ave))
        real_sens_mean.append(statistics.mean(c))"""
with open("sens_bias_total_random_cell_collective.txt","r") as fp:
    sens_bias=json.load(fp)
sens_bins=[[]]*15
sens_bias=[tuple(x) for x in sens_bias]
sens_bias=list(set(sens_bias))
for data in sens_bias:
    ind=int(np.floor((data[1]+0.035)*100/7))
    temp=sens_bins[ind].copy()
    temp.append(data[0])
    sens_bins[ind]=temp

sens_bias_yerr=[]

for bin in sens_bins:
    m=statistics.mean(bin)
    s=statistics.stdev(bin)
    sens_bias_yerr.append((m,s))
with open("sens_bias_total_random_cell_collective.txt","w") as fp:
    json.dump(sens_bias,fp)

with open("sens_bias_stdev_total_random_cell_collective.txt","w") as fp:
    json.dump(sens_bias_yerr,fp)
