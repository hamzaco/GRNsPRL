import pandas as pd
from Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures import Gene
import itertools
import os
import collections
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr

with open('dist_vs_input_number_brute_force.json', 'r') as fp:
    dist_vs_input_number_brute_force = json.load(fp)
with open('dist_vs_input_number_original.json', 'r') as fp:
    dist_vs_input_number_original = json.load(fp)
with open('dist_vs_input_number_symmetry_dictated.json', 'r') as fp:
    dist_vs_input_number_symmetry_dictated = json.load(fp)

for i in range(len(dist_vs_input_number_original)):
    lis=dist_vs_input_number_original[i]
    lis = [len(e) for e in lis]
    dist_vs_input_number_original[i]=lis

    lis = dist_vs_input_number_symmetry_dictated[i]
    lis = [len(e) for e in lis]
    dist_vs_input_number_symmetry_dictated[i] = lis

    lis = dist_vs_input_number_brute_force[i]
    lis = [len(e) for e in lis]
    dist_vs_input_number_brute_force[i] = lis

for input_number in range(2,12):

    freq = collections.Counter(dist_vs_input_number_original[input_number])
    freq=dict(freq)
    print(str(input_number)+' input:')
    normalization=sum(list(freq.values()))

    freq1= collections.Counter(dist_vs_input_number_symmetry_dictated[input_number])
    freq1 = dict(freq1)
    normalization1 = sum(list(freq1.values()))

    freq2 = collections.Counter(dist_vs_input_number_brute_force[input_number])
    freq2 = dict(freq2)
    normalization2 = sum(list(freq2.values()))
    # iterate dictionary named as freq to print
    # count of each element
    all_freq={**freq,**freq1,**freq2}
    all_freq=collections.OrderedDict(sorted(all_freq.items()))
    print(" {:^15} {:^7} {:^7} {:^7}".format('# of Classes', 'Original', 'SD', 'BF'))
    for key in all_freq:
        if key:
            if key not in freq:
                freq[key] = 0
            if key not in freq1:
                freq1[key] = 0
            if key not in freq2:
                freq2[key] = 0
            print(" {:<15} {:7.4f} {:7.4f} {:7.4f}".format(str(key),freq[key]/normalization, freq1[key]/normalization1 ,freq2[key]/normalization2))