from data_structures.Logic_Table_Network import Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network
from data_structures.Gene import  Gene
from scipy.stats import pearsonr
import itertools
from sample_functions import random_function



path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model"]

entries = os.listdir(path_cells)
all=[0]*20
not_ncf_o=[0]*20
asymmetric_o=[0]*20
not_canalizing_o=[0]*20
not_ncf_r=[0]*20
asymmetric_r=[0]*20
not_canalizing_r=[0]*20

for i in range(0,len(entries)):
    cell=entries[i]
    print(cell)
    if cell in big_entries:
        continue
    expressions = path + "/" + cell + "expressions.ALL.txt"
    externals = path + "/" + cell + "external_components.ALL.txt"
    original= Gene_Network(cell,expressions,externals)
    for gene in original.genes:
        if gene.number_of_inputs()<3 or gene.number_of_inputs()>8:
            continue
        nodes, is_ncf=gene.is_canalizing()
        all[gene.number_of_inputs()] += 1
        if not gene.is_symmetric():
            asymmetric_o[gene.number_of_inputs()] += 1
        if not is_ncf:
            not_ncf_o[gene.number_of_inputs()] += 1
        if len(nodes)==0:
            not_canalizing_o[gene.number_of_inputs()] += 1


input_numbers=[3,4,5,6,7,8]
sample_size=1000
for input_number in input_numbers:
    for i in range(sample_size):
        d=random_function(input_number)
        a=Gene()
        a.logical_dict=d
        a.in_degree=input_number
        a.is_influenced=set(map(lambda x: Gene(str(x)),np.arange(1,input_number+1)))
        nodes, is_ncf = a.is_canalizing()
        if not a.is_symmetric():
            asymmetric_r[a.number_of_inputs()] += 1
        if not is_ncf:
            not_ncf_r[a.number_of_inputs()] += 1
        if len(nodes)==0:
            not_canalizing_r[a.number_of_inputs()] += 1

print("input number \t\t\t\t 3 \t\t 4 \t\t 5 \t\t 6 \t\t 7 \t\t 8")
print("ncf percentages \t\t", str(round(1-not_ncf_o[3]/all[3],3)), "\t",
                            str(round(1-not_ncf_o[4] / all[4],3)), "\t",
                            str(round(1-not_ncf_o[5] / all[5],3)), "\t",
                            str(round(1-not_ncf_o[6] / all[6],3)), "\t",
                            str(round(1-not_ncf_o[7] / all[7],3)), "\t",
                            str(round(1-not_ncf_o[8] / all[8],3)), "\t",)

print("sym percentages \t\t", str(round(1-asymmetric_o[3]/all[3],3)), "\t",
                            str(round(1-asymmetric_o[4] / all[4],3)), "\t",
                            str(round(1-asymmetric_o[5] / all[5],3)), "\t",
                            str(round(1-asymmetric_o[6] / all[6],3)), "\t",
                            str(round(1-asymmetric_o[7] / all[7],3)), "\t",
                            str(round(1-asymmetric_o[8] / all[8],3)), "\t",)

print("can percentages \t\t", str(round(1-not_canalizing_o[3]/all[3],3)), "\t",
                            str(round(1-not_canalizing_o[4] / all[4],3)), "\t",
                            str(round(1-not_canalizing_o[5] / all[5],3)), "\t",
                            str(round(1-not_canalizing_o[6] / all[6],3)), "\t",
                            str(round(1-not_canalizing_o[7] / all[7],3)), "\t",
                            str(round(1-not_canalizing_o[8] / all[8],3)), "\t",)


print("ncf percentages random \t", str(round(1-not_ncf_r[3]/sample_size,3)), "\t",
                            str(round(1-not_ncf_r[4] / sample_size,3)), "\t",
                            str(round(1-not_ncf_r[5] / sample_size,3)), "\t",
                            str(round(1-not_ncf_r[6] / sample_size,3)), "\t",
                            str(round(1-not_ncf_r[7] / sample_size,3)), "\t",
                            str(round(1-not_ncf_r[8] / sample_size,3)), "\t",)

print("sym percentages random \t", str(round(1-asymmetric_r[3]/sample_size,3)), "\t",
                            str(round(1-asymmetric_r[4] / sample_size,3)), "\t",
                            str(round(1-asymmetric_r[5] / sample_size ,3)), "\t",
                            str(round(1-asymmetric_r[6] / sample_size,3)), "\t",
                            str(round(1-asymmetric_r[7] / sample_size,3)), "\t",
                            str(round(1-asymmetric_r[8] /sample_size,3)), "\t",)

print("can percentages random \t", str(round(1-not_canalizing_r[3]/sample_size,3)), "\t",
                            str(round(1-not_canalizing_r[4] /sample_size,3)), "\t",
                            str(round(1-not_canalizing_r[5] / sample_size,3)), "\t",
                            str(round(1-not_canalizing_r[6] / sample_size,3)), "\t",
                            str(round(1-not_canalizing_r[7] / sample_size,3)), "\t",
                            str(round(1-not_canalizing_r[8] / sample_size,3)), "\t",)

plt.figure()


percentages_canalizing_r=list(map(lambda x: 1-x/sample_size, not_canalizing_r[3:9]))
percentages_ncf_r=list(map(lambda x: 1-x/sample_size, not_ncf_r[3:9]))
percentages_symmetric_r=list(map(lambda x: 1-x/sample_size, asymmetric_r[3:9]))

percentages_canalizing_o=list(map(lambda x: 1-x, np.divide(not_canalizing_o[3:9],all[3:9])))
percentages_ncf_o=list(map(lambda x: 1-x, np.divide(not_ncf_o[3:9],all[3:9])))
percentages_symmetric_o=list(map(lambda x: 1-x, np.divide(asymmetric_o[3:9],all[3:9])))
plt.plot(np.arange(3,9),percentages_canalizing_r,label="canalizing")
plt.plot(np.arange(3,9),percentages_ncf_r,label="ncf")
plt.plot(np.arange(3,9),percentages_symmetric_r,label="symmetric")
plt.plot(np.arange(3,9),percentages_canalizing_o,label="canalizing original")
plt.plot(np.arange(3,9),percentages_ncf_o,label="ncf original ")
plt.plot(np.arange(3,9),percentages_symmetric_o,label="symmetric original")
plt.xlabel("Input number")
plt.ylabel("Percentage")
plt.legend()
plt.savefig("Frequency_of_features.png")
plt.show()