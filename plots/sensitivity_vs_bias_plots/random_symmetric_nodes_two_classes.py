import itertools
import random
import matplotlib.pyplot as plt
import statistics
import json
import os
import math
import pandas as pd
import altair as alt
import numpy as np
import collections
from data_structures.Logic_Table_Network import logical_gene,random_symmetric_gene
from scipy.stats import pearsonr
plt.figure()
max_input=10
sample_size=1000
count2n=0
got_it=0
type1={'Sensitivity':[],'Bias':[],'Distribution':[],'input_size':[],'small_class':[],'imbalance':[],'length_of_dist':[]}
for input_number in range(3,max_input):

    summ = 0
    distributions = []
    distributions_sens = []
    count1n = 0

    for i in range(sample_size):
        gene= random_symmetric_gene(input_number)
        s,real_size= gene.find_sensitivity()

        if real_size!= input_number:
            print(real_size)
            wwait=input('wtf: ')
        temp=0
        for value in list(gene.logical_dict.values()):
            if value:
                temp +=1
        bias=temp/(2**input_number)
        type1['Sensitivity'].append(s)
        type1['input_size'].append(real_size)
        type1['Bias'].append(bias)
        dist,_=gene.find_symmetries()
        dist= [int(e) for e in dist]
        minimum=min(dist)
        maximum=max(dist)
        imbalance=minimum/maximum
        dist.sort()
        dist= tuple(dist)
        type1['Distribution'].append(dist)
        type1['length_of_dist'].append(len(dist))
        type1['small_class'].append(minimum)
        type1["imbalance"].append(imbalance)
        if len(dist)==2:
            a = min(list(dist))
            guesses = []
            n=input_number
            count2n += 1

            """if  s==-2 * n / 2 ** n + 2 * a / 2 ** a + 2 * (n - a) / 2 ** (n - a) and\
                (bias==(2**a-1)*(2**(n-a)-1)/(2**n) or bias==1-(2**a-1)*(2**(n-a)-1)/(2**n)):
                got_it += 1
            elif s==2 * (n - 2 * a) / (2 ** n) + 2 * a / 2 ** a and\
                    (bias==(2**(n-a)-1)/(2**n) or bias==1-(2**(n-a)-1)/(2**n)):
                plt.plot(bias, s, '*k')
                got_it += 1
            elif s==-2 * (n - 2 * a) / (2 ** n) + 2 * (n - a) / 2 ** (n - a) and\
                    (bias==(2**(a)-1)/(2**n) or bias==1-(2**(a)-1)/(2**n)):
                plt.plot(bias, s, '*k')
                got_it += 1
"""


    """freq = collections.Counter(distributions)
    freq = dict(freq)
    print(str(input_number) + ' input:')
    normalization = sum(list(freq.values()))"""
plt.show()
#print((count1n+count2n)/(5*sample_size))
"""plt.xlabel('Activity bias')
plt.ylabel('Sensitivity')
plt.title(str(input_number)+'-inputs Genes')
    #plt.savefig(str(input_number)+'_input_bias_vs_sensitivity')
plt.show()"""
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"

with open(main_dir+'/dist_sens_bias_data/data_of_random_symmetric_xor_excluded_nodes.json', 'w') as fp:
    json.dump(type1, fp)

df=pd.DataFrame(type1)
alt.data_transformers.disable_max_rows()

chart=alt.Chart(df).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity',
    color='length_of_dist',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','imbalance','length_of_dist']
).interactive()

chart.show()