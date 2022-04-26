
from data_structures.Gene_Network import Gene_Network
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import math

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
entries = os.listdir(path_cells)

sens_activity_bias=[]
activity_bias=[]
c=0
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    for i in range(len(original.genes)):
        gene = original.genes[i]
        number_of_inputs= len(gene.is_influenced)+len(gene.is_influenced_externals)
        if number_of_inputs >13 or number_of_inputs==1:
            continue
        s,realsize=gene.find_sensitivity()
        if realsize==0:
            c += 1
            continue
        bias=gene.activity_bias()
        sens_activity_bias.append((s,bias))
print(c)
with open('sensitivity_data/sens_vs_bias.json', 'w') as fp:
    json.dump(sens_activity_bias, fp)

