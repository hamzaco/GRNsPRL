from data_structures.Logic_Table_Network import Logical_Table_Network
import json
import os
import numpy as np
import scipy.optimize as sp
from simulation import simulate
from data_structures.Gene_Network import Gene_Network, Gene
import statistics
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import random
import itertools

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",]
entries = os.listdir(path_cells)
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"

def random_partial_ncf(input_number,canalization_level,number_of_truths):
    truth=[True,False]
    reducible=True
    is_pcf=False
    a = Gene()
    if canalization_level>=input_number-1:
        return random_function(input_number,method="ncf")

    while reducible or not is_pcf:
        inputs = list(map(str, np.arange(1, input_number + 1)))
        canalizing_order = list(map(str, np.arange(1, canalization_level + 1)))
        random.shuffle(canalizing_order)
        lines = []
        not_2 = np.base_repr(number_of_truths, 2)
        if len(not_2) < canalization_level:
            not_2 = '0' * (canalization_level - len(not_2)) + not_2
        for i in range(len(canalizing_order)):
            canalizing_input = canalizing_order[i]
            lines.append((canalizing_input, random.choice(truth), not_2[i]))
        logical_dict = {}
        for c in itertools.product(truth, repeat=len(inputs)):
            flag=False
            for i in range(len(lines)):
                line = lines[i]

                canalizing_input = line[0]
                canalizing_value = line[1]
                canalized_value = line[2]
                if canalizing_value == c[i]:
                    logical_dict[c] = canalized_value=='1'
                    flag=True
                    break
            if not flag:
                logical_dict[c]=random.choice(truth)

        a.is_influenced=set([Gene(str(i)) for i in range(input_number)])
        a.logical_dict=logical_dict
        reducible= a.find_sensitivity()[1]!=input_number
        is_pcf= len(a.is_canalizing()[0])==canalization_level
    return logical_dict

def random_function(input_number,method=None,number_of_truths=0,sample_function=False):
    truth=[True,False]
    if method=="ncf":
        inputs=list(map(str,np.arange(1,input_number+1)))
        canalizing_order = inputs
        random.shuffle(canalizing_order)
        lines = []
        if number_of_truths>0:
            not_2=np.base_repr(number_of_truths,2)
            if len(not_2)<input_number:
                not_2= '0'*(input_number-len(not_2))+not_2
            for i in range(len(canalizing_order)):
                canalizing_input = canalizing_order[i]
                lines.append((canalizing_input, random.choice(truth), not_2[i]))  # input,canalized variable, canalized
            lines.append(("default", '0'))
        else:
            for i in range(len(canalizing_order) - 1):
                canalizing_input = canalizing_order[i]
                if sample_function:
                    lines.append((canalizing_input, random.choice(truth),
                                  random.random() < (np.exp(-2**(-i-1)*7)/(1+np.exp(-2**(-i-1)*7)) )))  # input,canalized variable, canalized value
                else:
                    lines.append((canalizing_input, random.choice(truth),
                                    random.choice(truth)))  # input,canalized variable, canalized value
            lines.append((canalizing_order[-1], random.choice(truth), True))
            lines.append(("default", False))

        logical_dict = {}
        for c in itertools.product(truth, repeat=len(canalizing_order)):
            for i in range(len(lines)):
                line = lines[i]

                canalizing_input = line[0]
                if canalizing_input == "default":
                    logical_dict[c] = line[1]
                    break
                canalizing_value = line[1]
                canalized_value = line[2]
                if canalizing_value == c[i]:
                    logical_dict[c] = canalized_value
                    break
        return logical_dict
    flag = True  # generate until condition is satisfied
    gene=Gene()
    gene.in_degree=input_number
    gene.is_influenced=set(map(lambda x: Gene(str(x)),np.arange(1,input_number+1)))
    while flag:
        count = 0
        logical_dict = {}
        if number_of_truths==0:
            number_of_truths = random.randint(1, (2 **input_number) - 1)
        truths = [True] * number_of_truths + [False] * (

                2 ** (input_number) - number_of_truths)
        random.shuffle(truths)
        for c in itertools.product(truth, repeat=input_number):
            logical_dict[c] = truths[count]
            count += 1
        gene.logical_dict = logical_dict
        s, real_size = gene.find_sensitivity()
        legit_methods = ["xor excluded", "symmetric", "asymmetric", "canalizing"]
        if not real_size == input_number:  # regenerate
            continue
        if method == "xor excluded" and (not gene.is_symmetric() or gene.behaviours()[0] > 0):
            continue
        if method == "symmetric" and (not gene.is_symmetric()):
            continue
        if method == "asymmetric" and gene.is_symmetric():
            continue
        if method == "canalizing" and len(gene.is_canalizing()[0]) == 0:
            continue
        if method == "nested canalizing" and not gene.is_canalizing()[1]:
            continue
        flag = False
    return logical_dict

def dict_to_table(dict):
    new_dict={}
    for key,value in dict.items():
        new_dict[tuple(map(int,key))]=int(value)
    return new_dict

sens=[]
bias=[]
sens_bias=[]
"""input_numbers=[2,3,4,5,6,7,8,9,10]
for input_number in input_numbers:
    sens = []
    for i in range(1,2**input_number):
        d=random_function(input_number,number_of_truths=i,method="ncf")
        a=Gene()
        a.logical_dict=d
        s,_=a.find_sensitivity()
        b= a.activity_bias()
        sens.append(s)
        bias.append(b)
        sens_bias.append((s,b))

with open("dist_sens_bias_data/sens_random_ncf_equal_prob.txt","w") as fp:
    json.dump(sens, fp)

plt.hist(bias,bins=20)
plt.show()"""
