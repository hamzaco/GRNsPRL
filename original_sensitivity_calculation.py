import json
import os

import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
import itertools

sens_of_networks={}
path = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"


entries = os.listdir(path)
dist_sens_bias=[[]]*14
sens_bias=[]
for i in range(0,len(entries)):
    cell=entries[i]
    temp=[]
    print(cell + ":")
    expressions = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"+ "/" + cell + "expressions.ALL.txt"
    externals = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered" + "/" + cell + "external_components.ALL.txt"
    gn = Gene_Network(cell, expressions, externals)
    for gene in gn.genes:
        if gene.number_of_inputs()>5 or gene.number_of_inputs()==1:
            continue
        s,real_size=gene.find_sensitivity()
        if s > 0:
            bias=gene.activity_bias()
            dist,_=gene.find_symmetries()
            temp=dist_sens_bias[real_size].copy()
            temp.append((dist,s,bias))
            sens_bias.append((s,bias))
            dist_sens_bias[real_size]=temp


"""with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original.json','w') as fp:
    json.dump(dist_sens_bias,fp)
"""
with open(main_dir+'/sensitivity_data/sens_vs_bias_upto4.json','w') as fp:
    json.dump(sens_bias,fp)