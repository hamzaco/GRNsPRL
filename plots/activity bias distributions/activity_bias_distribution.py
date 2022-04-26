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

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",]
entries = os.listdir(path_cells)
biases=[[]]*20
for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for gene in original.genes:
        t=biases[gene.number_of_inputs()].copy()
        t.append(gene.activity_bias())
        biases[gene.number_of_inputs()]=t

yerr=[]
y=[]
for input_number in range(1,12):
    data=biases[input_number]
    yerr.append(statistics.stdev(data))
    y.append(statistics.mean(data))


