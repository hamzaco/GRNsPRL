
from data_structures.Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import statistics
import random
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from data_structures.Gene import is_symmetric
import itertools
import math
def probability_of_symmetry(dist):
    total=sum(dist)
    indexs=[i for i in range(total)]
    partitions=[]
    for p in dist:
        temp=indexs[0:p]
        indexs=indexs[p:]
        partitions.append(temp)
    indexs=[i for i in range(total)]
    iteration_size=1
    symmetries=0
    for i in range(iteration_size):
        random.shuffle(indexs)
        pair=indexs[0:2]
        for p in partitions:
            if pair[0] in p and pair[1] in p:
                symmetries +=1
    return symmetries/iteration_size

def probability_of_consistency(dists):
    p=1
    ps=1
    for dist in dists:
        temp= probability_of_symmetry(dist)
        p *= temp
        ps *= (1-temp)
    return p+ps


def intersection(lst1, lst2):
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

def get_expectation_consistency(genes,dist_dict):
    expected=0
    expected_sampled=0
    var_sampled=0
    all=0
    in_reality=0
    for pair in itertools.combinations(genes, 2):
        Cij=0
        s = 1
        not_s = 1
        symmetry_list=[]
        dists=[]
        for gene in genes:
            inputs = list(gene.is_influenced) + list(gene.is_influenced_externals)
            n=len(inputs)
            if n<2:
                continue
            if pair[0] in inputs and pair[1] in inputs:
                Dij = 0
                dist=dist_dict[gene]
                dists.append(dist)
                for p in dist:
                    if p>1:
                        Dij += nCr(p,2)/nCr(n,2)
                s *= Dij
                not_s *= 1-Dij
                symmetry_list.append(is_symmetric(pair[0].name, pair[1].name,gene)[0])

        Cij= s+not_s
        samples=[]
        sample_size=100
        Cij_sampled = probability_of_consistency(dists)
        if Cij_sampled != 2:
            for i in range(sample_size):
                samples.append(probability_of_consistency(dists))
            var_sampled += statistics.variance(samples)
            expected_sampled += statistics.mean(samples)

        if Cij!=2:
            all +=1
            expected += Cij
            in_reality += symmetry_list==[True]*len(symmetry_list) or symmetry_list==[False]*len(symmetry_list)

    return expected,expected_sampled,var_sampled,in_reality, all
def nPr(n,r):
    return math.factorial(n)/(math.factorial(n-r))
def nCr(n,r):
    return math.factorial(n)/(math.factorial(n-r)*math.factorial(r))
def NumberOfInputs(Genes):
    l = []
    for gene in Genes:
       l.append((len(gene.is_influenced),len(gene.is_influenced_externals)))
    return l
path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]

entries = os.listdir(path_cells)
symmetric=[]
cross_symmetric=[]
non_symmetric=[]
temp=0
temp1=0
total_expected_inconsistent=0
all_pairs=0
total_inconsistent=0
for i in range(0,len(entries)):
    cell=entries[i]
    if cell in big_entries:
        continue
    print(cell + ":")
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    inputs = NumberOfInputs(original.genes)
    symmetry_private = []
    non_symmetry_private = []
    dist_dict={}
    for gene in original.genes:
        dist,_=gene.find_symmetries()
        dist_dict[gene]=dist

    expected,expected_s,var_s,in_reality,all_related_pairs=get_expectation_consistency(original.genes,dist_dict)
    print(expected,expected_s,var_s,all_related_pairs,in_reality)
    all_pairs += all_related_pairs
    total_expected_inconsistent += all_related_pairs-expected

    total_inconsistent += all_related_pairs-in_reality
    """for gene in original.genes:
        inputs_genes = [i.name for i in list(gene.is_influenced)]
        inputs = inputs_genes+ list(gene.is_influenced_externals)

        n= len(inputs)

        if n<2:
            continue

        for pair in itertools.combinations(inputs_genes,2):
            is_sym,_,cross_sym=is_symmetric(pair[0],pair[1],gene)
            if  is_sym:
                symmetric.append(pair)
                symmetry_private.append(pair)
                if cross_sym:
                    cross_symmetric.append(pair)
            else:
                non_symmetric.append(pair)
                non_symmetry_private.append(pair)

    symmetry_private = set(symmetry_private)
    non_symmetry_private = set(non_symmetry_private)
    intersect=symmetry_private.intersection(non_symmetry_private)

    print(len(symmetry_private),len(non_symmetry_private),len(intersect))
    if len(symmetry_private)>0:

        temp += len(intersect)/len(set(symmetry_private))
    if len(non_symmetry_private)>0:
        temp1 += len(intersect)/len(set(non_symmetry_private))"""
"""with open('symmetic_relations.json', 'w') as fp:
    json.dump(symmetric,fp)

with open('non_symmetic_relations.json', 'w') as fp:
    json.dump(non_symmetric,fp)
with open('cross_symmetic_relations.json', 'w') as fp:
    json.dump(cross_symmetric,fp)
"""

with open('general_data/symmetic_relations.json', 'r') as fp:
    symmetric = json.load(fp)

with open('general_data/non_symmetic_relations.json', 'r') as fp:
    non_symmetric = json.load(fp)


with open('general_data/cross_symmetic_relations.json', 'r') as fp:
    cross_symmetric = json.load(fp)
symmetric=set([tuple(sorted(a)) for a in symmetric])
non_symmetric=set([tuple(sorted(a))  for a in non_symmetric])
cross_symmetric=set([tuple(sorted(a)) for a in cross_symmetric])

intersect=non_symmetric.intersection(symmetric)
print("Symmetric: ",len(symmetric),"Non-symmetric: ",len(non_symmetric),"Inconsistnet: ",len(intersect))
print("Symmetric pairs: ",len(symmetric),"Cross-symmetric: ",len(cross_symmetric), "percentage: ",
      len(cross_symmetric)/len(symmetric))

print("Expected number of inconsistent pairs:", total_expected_inconsistent)
print("Total number of inconsistent pairs:", total_inconsistent)

print("All related pairs: ", all_pairs)
