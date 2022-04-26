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
input_numbers=[2,3,4,5,6]
d=[0,0]
path="C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
for input_number in  input_numbers:
    truths=[0]*(2**(input_number))
    for i in range(0,len(entries)):
        cell=entries[i]
        if cell in big_entries:
            continue
        expressions = path + "/" + cell + "expressions.ALL.txt"
        externals = path + "/" + cell + "external_components.ALL.txt"
        gn = Gene_Network(cell, expressions, externals)
        #gn.sensitivity_test()
        for gene in gn.genes:
            if len(gene.is_influenced)+len(gene.is_influenced_externals)!=input_number:
                continue
           # temp,classes,R=gene.find_symmetries()
            truth=gene.activity_bias()*(2**input_number)
            truths[int(truth)] +=1

           # if temp in distributions_and_sensitivities:
            #    a= list(distributions_and_sensitivities[temp])
             #   a.append(gene.sensitivity)
              #  a=tuple(a)
               # distributions_and_sensitivities[temp]=a
            #else:
             #   distributions_and_sensitivities[temp]=(gene.sensitivity,gene.sensitivity)

            if truth in truths_and_sensitivities:
                a= list(truths_and_sensitivities[truth])
                a.append(gene.sensitivity)
                a=tuple(a)
                truths_and_sensitivities[truth]=a
            else:
                truths_and_sensitivities[truth]=(gene.sensitivity,gene.sensitivity)
    plt.figure()
    i=0
    #for distribution in distributions_and_sensitivities:
     #   for point in distributions_and_sensitivities[distribution]:
      #      plt.plot(distribution,point,'*')
       # i +=1

    plt.figure()
    i=0
    #for distribution in truths_and_sensitivities:
     #   for point in truths_and_sensitivities[distribution]:
      #      plt.plot(distribution,point,'*k')
       # i +=1
    plt.show()

    plt.figure()
    for i in range(len(truths)):
        plt.plot(i,truths[i],'*k')
        if i%2==0:
            d[0]+=truths[i]
        else:
            d[1] += truths[i]
    plt.grid()
    plt.xlabel('Number of Truths')
    plt.ylabel('Number of genes')
    plt.title('Number of Truths in Full Logical Table for '+str(input_number)+'-input')
    plt.show()


plt.bar(['even','odd'], d, width=0.8, align='center', )
plt.title('Number of Truths for all functions')
plt.savefig('number_of_truths_all')
plt.show()

#with open('distributions_and_sensitivities.json', 'w') as fp:
 #   json.dump(distributions_and_sensitivities, fp)
#with open('truths_and_sensitivities.json', 'w') as fp:
 #   json.dump(truths_and_sensitivities, fp)