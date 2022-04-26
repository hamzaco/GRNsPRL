
from data_structures.Logic_Table_Network import Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network
from scipy.stats import pearsonr
import itertools

path = "C:/Users/Lenovo/PycharmProjects/INDEP"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
entries = os.listdir(path_cells)
big_entries= ["CC_EGFR_ErbB_Signaling",
              "CC_CD4p_TCell_Diff_Plas",
              "CC_HIV1_interactions_with_TCell_Signalling_Pathway",
"CC_Inflammatory_Bowel_Disease_IBD_Model"]

def sensitivities_random_network(method):
    with open('sensitivity_data/sensitivities.json', 'r') as fp:
        sensitivities = json.load(fp)

    partitions=[[]]*20
    activity_bias=[]
    random_sensitivities={}

    for i in range(0,len(entries)):
        cell=entries[i]
        if cell in big_entries:
            continue
        print(cell + ":")
        expressions = path + "/network_data_altered/" + cell + "expressions.ALL.txt"
        externals = path + "/network_data_altered/" + cell + "external_components.ALL.txt"
        original= Gene_Network(cell,expressions,externals)
        gn = Logical_Table_Network(cell, expressions, externals)
        gn.readData()
        random_size=1
        sensitivities_data=[]
        for i in range(0,random_size):

            rand2 = random_network(original,method)
            summ = 0
            for gene in rand2.genes:
                dist,_= gene.find_symmetries()
                s, real_size = gene.find_sensitivity()
                if real_size != len(gene.is_influenced)+len(gene.is_influenced_externals):
                    print(gene)
                    w=input("waiting:  ")
                summ += s
                bias=gene.activity_bias()
                temp=partitions[real_size].copy()
                temp.append((dist,s,bias))
                partitions[real_size]=temp
                activity_bias.append(bias)
            sensitivities_data.append(summ/len(rand2.genes))
        random_sensitivities[cell]=sum(sensitivities_data)/len(sensitivities_data)
        print(random_sensitivities[cell])

    with open('dist_sens_bias_data/dist_sens_bias_'+method+'1.json', 'w') as fp:
        json.dump(partitions, fp)

    sensitivities=list(sensitivities.values())
    random_sensitivities= list(random_sensitivities.values())

    #plt.figure()
    """plt.hist(sensitivities,bins=30,range=[0, 3],weights=[1/len(sensitivities)]*len(sensitivities),
             facecolor='red', ec='black', alpha=0.2, label='Original')"""
    """plt.hist(random_sensitivities,bins=30,range=[0, 3],
             weights=[1/len(random_sensitivities)]*len(random_sensitivities),
             facecolor='blue', ec='black', alpha=0.2, label='Random NCF')"""
    """plt.hist(symmetric_sensitivities,bins=30,range=[0, 3],
             weights=[1/len(symmetric_sensitivities)]*len(symmetric_sensitivities),
             facecolor='green', ec='black', alpha=0.2, label='Symmetry dictation')
    
    plt.hist(asymmetric_sensitivities,bins=30,range=[0, 3],
             weights=[1/len(asymmetric_sensitivities)]*len(asymmetric_sensitivities),
             facecolor='orange', ec='black', alpha=0.2, label='Asymmetric')"""

    plt.legend()
    #plt.savefig('Sensitivities Histogram Random without xor')
    plt.show()


    """
    n_brute_force,bins, _= plt.hist(activity_bias_original, bins=20, range=[0, 1],
            facecolor='red',ec='black', alpha=0.2,
            weights=[1/len(activity_bias_original)]*len(activity_bias_original),
                                    label='Original')
    n_symmetry_dictated,_,_=plt.hist(activity_bias_random_symmetry_dictated, bins=20,
            range=[0, 1],facecolor='blue', ec='black',alpha=0.2,
        weights=[1/len(activity_bias_random_symmetry_dictated)]*len(activity_bias_random_symmetry_dictated),
                                     label='Symmetry dictated')
    n_original,_,_=plt.hist(activity_bias_random_brute_force, bins=20, range=[0, 1],
            facecolor='green', ec='black',alpha=0.2,
            weights=[1/len(activity_bias_random_brute_force)]*len(activity_bias_random_brute_force) ,
                            label='Brute force')
    plt.xlabel('Activity bias')
    plt.legend()
    plt.savefig('Activity bias distributions with two random method')
    plt.show()
    n_asymmetric,_,_=plt.hist(activity_bias_random_asymmetric, bins=20, range=[0, 1],
            facecolor='green', ec='black',alpha=0.2,
            weights=[1/len(activity_bias_random_asymmetric)]*len(activity_bias_random_asymmetric) ,
                            label='Asymmetric')
    
    
    corr, _ = pearsonr(n_symmetry_dictated,n_original)
    print(corr)
    corr, _ = pearsonr(n_brute_force,n_original)
    print(corr)
    corr, _ = pearsonr(n_asymmetric,n_original)
    print(corr)"""



sensitivities_random_network(method="ncf")