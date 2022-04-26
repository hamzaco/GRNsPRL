
from data_structures.Logic_Table_Network import compare_logic_dict, restrict
import json
import os
import math
import numpy as np
import statistics
import matplotlib.pyplot as plt
import scipy.stats as st
import random
import altair as alt
import pandas as pd
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from data_structures.Gene import is_symmetric
import itertools

def shuffle_in_degrees(Network):
    for g in Network.genes:
        g.shuffle_in_degree(Network.genes,Network.externals)

def get_consistency(network): #genes is a list of genes for a particular network.
    genes=network.genes
    externals=network.externals
    consistent=0
    not_consistent=0
    all_inputs=[g.name for g in genes] # +list(externals)
    for pair in itertools.combinations(all_inputs, 2): #iterate over all pairs
        symmetry_list=[] #store the information about the symmetry
        for gene in genes: #iterate over all possible genes in the network
            inputs = [g.name for g in gene.is_influenced]+ list(gene.is_influenced_externals)
            if pair[0] in inputs and pair[1] in inputs:
                symmetry_list.append(is_symmetric(pair[0], pair[1],gene)[0]) #check whether they are in the same symmetry class
        if len(symmetry_list)>1:
            if symmetry_list==[True]*len(symmetry_list) or symmetry_list==[False]*len(symmetry_list): #check consistency, they are always symmetric or asymmetric.
                consistent += 1
            else:
                not_consistent += 1
    return consistent,not_consistent


def nPr(n,r):
    return math.factorial(n)/(math.factorial(n-r))
def nCr(n,r):
    return math.factorial(n)/(math.factorial(n-r)*math.factorial(r))
path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]
networks=[]
entries = os.listdir(path_cells)
for i in range(0, len(entries)):
    cell = entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    networks.append(original)
    for gene in original.genes:
        _,classes=gene.find_symmetries()
        gene.classes=classes

network_samples=[]
iteration_size=100
for j in range(iteration_size):
    sum_c = 0
    sum_ic = 0
    for i in range(0,len(networks)):
        original= networks[i]
        if j>0:
            shuffle_in_degrees(original)
        c,ic=get_consistency(original)
        sum_c += c
        sum_ic += ic
    print(sum_c,sum_ic)
    network_samples.append(sum_ic/(sum_ic+sum_c))
print(statistics.mean(network_samples))
print(statistics.stdev(network_samples))