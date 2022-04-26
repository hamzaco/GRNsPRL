from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from Gene_Network import Gene_Network, shuffle_network,print_genes
import itertools
def is_same_cycle(eq, cycle):
    for i in range(0,len(cycle)):
        temp = cycle[1:]
        temp.append(cycle[0])
        cycle=temp
        if eq == cycle:
            return True
    return False
path = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
big_entries= ["CC_CD4p_Tcell_Diff",
"CC_Differentiation_T_lymphocytes",
"CC_Inflammatory_Bowel_Disease_IBD_Model",
"CC_Glucose_Repression_Signaling_2009",
"CC_HCC1954_Breast_Cell_Line_Long-term_ErbB_Network",
"CC_IGVH_mutations",
"CC_Lymphopoiesis_Regulatory_Network",
"CC_Signaling_in_Macrophage_Activation",
"CC_Signal_Transduction_in_Fibroblasts",
"CC_Treatment_of_Castration_Resistant_Prostate_Cancer",
"CC_Yeast_Apoptosis"]
all_attractors={}
corresponding_ext={}
entries = os.listdir(path)
transition_distribution={}
with open('corresponding_externals.json', 'r') as fp:
    corresponding_ext = json.load(fp)
with open('attractors.json', 'r') as fp:
    all_attractors = json.load(fp)

entries = os.listdir(path)
sens_over_fixed_point={}
transition_distribution = {}
for i in range(0, len(entries)):

    cell = entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "/" + "expressions.ALL.txt"
    externals = path + "/" + cell + "/" + "external_components.ALL.txt"
    grn = Gene_Network(cell, expressions, externals, attractors=all_attractors[cell],
                       corresponding_externals=corresponding_ext)
    grn.initialize_effects()
    sens = 0
    weight=0
    corresponding_externals=corresponding_ext[cell]
    for attractor in grn.attractors:
        if isinstance(attractor[0], str):
            for i in range(0, len(attractor[0])):
                if attractor[0][i] == '1':
                    g, conditions = grn.set_conditions(attractor[0][:i] + '0' + attractor[0][i + 1:],
                                                       exstate=corresponding_externals[0])
                else:
                    g, conditions = grn.set_conditions(attractor[0][:i] + '1' + attractor[0][i + 1:],
                                                       exstate=corresponding_externals[0])
                eq = grn.find_equilibrium(conditions)
                if eq != attractor[0]:
                    sens += 1
                weight +=1
        else:
            for i in range(0, len(attractor[0][0])):
                if attractor[0][0][i] == '1':
                    g, conditions = grn.set_conditions(attractor[0][0][:i] + '0' + attractor[0][0][i + 1:],
                                                       exstate=corresponding_externals[0])
                else:
                    g, conditions = grn.set_conditions(attractor[0][0][:i] + '1' + attractor[0][0][i + 1:],
                                                       exstate=corresponding_externals[0])

                eq = grn.find_equilibrium(conditions)
                if isinstance(eq, str):
                    sens += 1
                elif is_same_cycle(eq, attractor[0]):
                    sens += 1
                weight += 1
        corresponding_externals.remove(corresponding_externals[0])
    print(sens/weight)
    sens_over_fixed_point[cell]=sens/weight

with open('sens_over_fixed_point.json', 'w') as fp:
    json.dump(sens_over_fixed_point, fp)
