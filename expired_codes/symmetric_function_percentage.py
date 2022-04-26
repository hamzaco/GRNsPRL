
import itertools
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
import statistics
import json
import os
from data_structures.Gene import Gene
import math
from Logic_Table_Network import logical_gene
from scipy.stats import pearsonr
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
truth=[True,False]
max_input=10
iteration_size=1000
plt.figure()
for input in range(2,max_input):
    odd_numbers=[]
    print(input)
    symmetric_count=0
    effectivity_count=0
    for i in range(1,2**(input)):
        odd_numbers.append(i)
    count = iteration_size
    for j in range (0,iteration_size):
        logical_dict={}
        number_of_truths=random.choice(odd_numbers)
        truths = [True] * number_of_truths + [False] * (2 **input - number_of_truths)

        random.shuffle(truths)
        i=0
        for c in itertools.product(truth, repeat=input):
            logical_dict[c]=truths[i]
            i+=1
        gene=Gene(name='random '+str(j), function=logical_dict)

        dist, classes = gene.find_symmetries()
        if dist != list('1' * input):
            symmetric_count +=1
        s, real_size=gene.find_sensitivity()
        if real_size != input:
            effectivity_count +=1
    #plt.plot(input,symmetric_count/count,'*k')
    p=1-(1-0.5**(2**(input-2)))**nCr(input,2)
    #plt.plot(input,p,'*r')
    effectivity=(1-0.5**(2**(input-1)))**input
    plt.plot(input,effectivity_count/count,'*g')
    plt.plot(input,1-effectivity,'*c')

    symmetry_count = iteration_size - symmetric_count
    p = 2**(3*2**(input - 2))
    print(iteration_size*symmetry_count/p)

    #error_sym=  (symmetric_count/count-p)/(symmetric_count/count)

plt.show()