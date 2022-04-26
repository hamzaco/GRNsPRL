from Gene_Network import Gene_Network
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np

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
entries = os.listdir(path)

distributions_and_sensitivities={}
truths_and_sensitivities={}
input_numbers=[3]
d=[0,0]
for i in range(0,len(entries)):
        cell=entries[i]
        if cell not in big_entries:
            continue
        expressions = path + "/" + cell + "/" + "expressions.ALL.txt"
        externals = path + "/" + cell + "/" + "external_components.ALL.txt"
        gn = Gene_Network(cell, expressions, externals)
        gn.sensitivity_test()
        expressions_gene=open(cell+'expressions.ALL.txt','w')
        for gene in gn.genes:
            expressions_gene.write(gene.name+' = '+gene.function+'\n')
        expressions_gene.close()

        externals_gene=open(cell+'external_components.ALL.txt','w')
        for ext in gn.externals:
            externals_gene.write(ext+'\n')

        externals_gene.close()