from data_structures.Logic_Table_Network import Logical_Table_Network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, Gene
from scipy.stats import pearsonr
import random
import itertools
from sample_functions import random_function


input_numbers=[12]
sample_size=1000
biases=[]
sens=[]

for i in input_numbers:
    for j in range(1,2**i):
        a=Gene()
        a.logical_dict=random_function(i,method="ncf",number_of_truths=j)
        if (a.activity_bias()*2**i)%2==0:
            sens.append(a.find_sensitivity()[0])
            biases.append(a.activity_bias())


fig, axs = plt.subplots(2, 2)
axs[1,1].hist(biases,bins=20,ec="black",range=[0, 1])

axs[0,0].hist(sens,bins=20,ec="black",range=[0,1.25],orientation="horizontal")

axs[0,1].plot(biases,sens,"*")
plt.show()
print(sum(sens)/len(sens))

ind=np.argsort(biases)

biases_ncf=np.take_along_axis(np.array(biases),ind,axis=0)
sensitivities_ncf=np.take_along_axis(np.array(sens),ind,axis=0)

i=0
s=0
while s<1:
    i +=1
    s=sensitivities_ncf[i]
print(biases_ncf[i],s)

main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
with open(main_dir+'/sensitivity_data/sens_vs_bias.json', 'r') as fp:
    sensitivity_vs_acitivity_bias = json.load(fp)

"""sens_original=[]
bias_original=[]
for data in sensitivity_vs_acitivity_bias:
    sens_original.append(data[0])
    bias_original.append(data[1])

plt.hist(bias_original,bins=20,ec="black")
plt.show()

with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original.json', 'r') as fp:
    data_original = json.load(fp)

for data in data_original:
    if sum(data[0])>2:
        sens_original.append(data[1])
        bias_original.append(data[2])

plt.hist(bias_original,bins=20,ec="black")
plt.show()"""
