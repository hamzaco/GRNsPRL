from __future__ import print_function
from itertools import product       # forms cartesian products
import statistics
from matplotlib import pyplot as plt
from data_structures.Gene import Gene
import math
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
max_input=6
for n in range(2,max_input):
    asymmetry_count=0
    effectivity_count=0
    sens_low=[]
    sens_high=[]
    count=(2**(2**n))
    print('All possible truth tables for n =', n)
    inputs = list(product([False, True], repeat=n))
    for output in product([False, True], repeat=len(inputs)):
        logical_dict={}
        number_of_truths=0
        for row, result in zip(inputs, output):
            logical_dict[tuple(row)]=result
            if result:
                number_of_truths +=1
        gene=Gene('random',function=logical_dict)
        s, real_size=gene.find_sensitivity()
        dist,_=gene.find_symmetries()
        bias=number_of_truths/(2**n)
        if real_size != n:
            effectivity_count +=1

        if dist==[1]*n:
            asymmetry_count +=1
            """print(gene.logical_dict)
            wait=input('waiting:')"""
    print(count-asymmetry_count)
    """effectivity = (1 - 0.5 ** (2 ** (n - 1))) ** n
    print(1-effectivity)
    print(effectivity_count/count)"""

    """print(statistics.mean(sens_low))
    print(statistics.mean(sens_high))
    print(sens_low)
    print(sens_high)
    plt.figure()
    plt.hist(sens_low,bins=40, range=[0,4])
    plt.show()"""

