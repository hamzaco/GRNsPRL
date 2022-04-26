from data_structures.Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from pylab import cm
import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

with open(main_dir+'/sensitivity_data/sens_vs_bias.json', 'r') as fp:
    sensitivity_vs_acitivity_bias = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_ncf.json', 'r') as fp:
    data_ncf = json.load(fp)

coverages=[]
deviations=np.arange(0,1.01,step=0.005)

for deviation in deviations:
    c=0
    in_range=0
    for key in sensitivity_vs_acitivity_bias:
        s=key[0]
        c += 1
        #c += freq_table[key]
        if abs(s-1)<=deviation:
            in_range += 1
            #in_range += freq_table[key]

    coverages.append(in_range/c)
plt.figure()
plt.plot(deviations,coverages,label="Cell Collective")

coverages=[]
deviations=np.arange(0,1.01,step=0.005)

with open(main_dir+"/dist_sens_bias_data/sens_random_ncf_equal_prob.txt", "r") as fp:
    sens_ncf_random=json.load(fp)
for deviation in deviations:
    c=0
    in_range=0
    for key in sens_ncf_random:
        s=key
        c += 1
        #c += freq_table[key]
        if abs(s-1)<deviation:
            in_range += 1
            #in_range += freq_table[key]

    coverages.append(in_range/c)
plt.plot(deviations,coverages,label="Random NCF",c="k")


coverages=[]
deviations=np.arange(0,1.01,step=0.005)
with open(main_dir+"/dist_sens_bias_data/sens_random_3.txt", "r") as fp:
    sens_ncf_random=json.load(fp)
for deviation in deviations:
    c=0
    in_range=0
    for key in sens_ncf_random:
        s=key
        c += 1
        #c += freq_table[key]
        if abs(s-1)<deviation:
            in_range += 1
            #in_range += freq_table[key]
    coverages.append(in_range/c)
#plt.plot(deviations,coverages,c="purple")


coverages=[]
deviations=np.arange(0,1.01,step=0.005)
with open(main_dir+"/dist_sens_bias_data/sens_random_4.txt", "r") as fp:
    sens_ncf_random=json.load(fp)
for deviation in deviations:
    c=0
    in_range=0
    for key in sens_ncf_random:
        s=key
        c += 1
        #c += freq_table[key]
        if abs(s-1)<deviation:
            in_range += 1
            #in_range += freq_table[key]
    coverages.append(in_range/c)


#plt.plot(deviations,coverages,c="g")
coverages=[]
deviations=np.arange(0,1.01,step=0.005)
with open(main_dir+"/dist_sens_bias_data/sens_random_6.txt", "r") as fp:
    sens_ncf_random=json.load(fp)
for deviation in deviations:
    c=0
    in_range=0
    for key in sens_ncf_random:
        s=key
        c += 1
        #c += freq_table[key]
        if abs(s-1)<deviation:
            in_range += 1
            #in_range += freq_table[key]
    coverages.append(in_range/c)
#plt.plot(deviations,coverages,c="orange")

coverages=[]
deviations=np.arange(0,1.01,step=0.05)
with open(main_dir+"/dist_sens_bias_data/sens_random_10.txt", "r") as fp:
    sens_ncf_random=json.load(fp)
for deviation in deviations:
    c=0
    in_range=0
    for key in sens_ncf_random:
        s=key
        c += 1
        #c += freq_table[key]
        if abs(s-1)<deviation:
            in_range += 1
            #in_range += freq_table[key]

    coverages.append(in_range/c)
#plt.plot(deviations,coverages,c="r")

plt.xlim([0,1])
plt.ylim([0,1])
#plt.text(0.85,0.88,"n=3",fontsize=15)
#plt.text(0.85,0.68,"n=4",fontsize=15)
#plt.text(0.85,0.4,"n=6",fontsize=15)
#plt.text(0.85,0.22,"n=10",fontsize=15)
plt.xticks(np.arange(0,1.005,step=0.2),fontsize=15)
plt.yticks(np.arange(0,1.005,step=0.2),fontsize=15)
plt.grid()
plt.xlabel(r'Distance to $\xi=1$ ',fontsize=20)
plt.ylabel("Probability",fontsize=20)

#plt.legend(fontsize=9)
plt.legend(loc=(0.71,0))
plt.title("\n",fontsize=20)
#plt.text(0.45,1.07,"(b)",fontsize=20)
plt.tight_layout()
plt.savefig("cov_pre_random_ncf_equal_probability_reduced")
plt.show()