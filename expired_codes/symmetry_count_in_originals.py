from Gene_Network import Gene_Network
from Logic_Table_Network import Logical_Table_Network, logical_gene,random_network
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import math

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
path = "C:/Users/Lenovo/PycharmProjects/INDEP"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",
]
entries = os.listdir(path_cells)


input_numbers=[3]
d=[0,0]

differences_in_symmetry={}

count=[0]*20
count_random=[0]*20
all=[0]*20
for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original = Gene_Network(cell, expressions, externals)
    symmetry_counts = [0]*len(original.genes)
    for i in range(0,len( original.genes)):
        gene= original.genes[i]
        dist= gene.find_symmetries()
        number_of_inputs=len(gene.is_influenced_externals)+len(gene.is_influenced)
        all[number_of_inputs] +=1
        if number_of_inputs<=2:
            continue
        symmetry_count=0
        for class_size in dist:

            if class_size > 1:
                symmetry_count += nCr(class_size,2)
        if dist==int('1'*number_of_inputs):
            count[number_of_inputs]=+1
        symmetry_counts[i]=symmetry_count/nCr(number_of_inputs,2)
    random_size=10
    gn = Logical_Table_Network(cell, expressions, externals)
    gn.readData()
    symmetry_counts_random = [0]*len(symmetry_counts)

    for j in range(0, random_size):
        gn = random_network(gn, method="brute force")
        for i in range(0,len(gn.genes)):
            gene= gn.genes[i]
            ##symmetry_count=0
            #dist, classes = gene.find_symmetries()
            number_of_inputs = len(gene.is_influenced_externals) + len(gene.is_influenced)
            all[number_of_inputs] += 1
            if number_of_inputs <= 2:
                continue
            for class_size in dist:
                if class_size > 1:
                    symmetry_count += nCr(class_size, 2)
            if dist == int('1' * number_of_inputs):
                count_random[number_of_inputs] = +1
            symmetry_counts_random[i] += symmetry_count/random_size/nCr(number_of_inputs,2)
    difference=[0]*len(symmetry_counts)
    for i in range(0,len(symmetry_counts)):
        difference[i]=(symmetry_counts[i]-symmetry_counts_random[i])
    differences_in_symmetry[cell]=difference
    print(difference)
    print(symmetry_counts)
    print(symmetry_counts_random)
with open('differences_in_symmetry_totally random.json', 'w') as fp:
    json.dump(differences_in_symmetry, fp)