from data_structures.Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures.Gene import Gene
import itertools
from matplotlib import pyplot as plt
import os
from data_structures.BooleanNetwork import *
from sample_functions import random_partial_ncf
sample_size=1000

for inputs in range(7,16):
    sens = []
    bias = []
    if inputs==8:
        continue
    odds = np.arange(1, 2 **7, 1)
    print(inputs)
    for i in range(sample_size):
        a=Gene()
        a.logical_dict=random_partial_ncf(inputs,7,random.choice(odds))
        sens.append(a.find_sensitivity()[0])
        bias.append(a.activity_bias())

    plt.plot(bias,sens,"*",markersize=0.2,label=str(inputs))
plt.legend()
plt.savefig("pcf.pdf")
plt.show()

def beta2(n,alpha):
    if n==4:
        return random.randint(12,48)
    if alpha < 2**(n-1):
        return beta2(n-1,alpha)+alpha
    return beta2(n-1,2**n-alpha)+(2**n-alpha)

def beta(n,alpha):
    if n==1 and alpha==1:
        return 1
    if alpha < 2**(n-1):
        return beta(n-1,alpha)+alpha
    return beta(n-1,2**n-alpha)+(2**n-alpha)


"""n=9
evens=np.arange(1,2**n,1)
for e in evens:
    sens.append(beta2(n,e)/(2**(n-1)))
    bias.append(e/(2**n))

plt.plot(bias,sens,"o",label="n=9", markersize=3)

sens=[]
bias=[]

n=12
evens=np.arange(1,2**n,1)
for e in evens:
    sens.append(beta2(n,e)/(2**(n-1)))
    bias.append(e/(2**n))

plt.plot(bias,sens,"o",label="n=12", markersize=0.2)


sens=[]
bias=[]
odds=np.arange(1,2**n,2)
for e in odds:
    sens.append(beta(n,e)/(2**(n-1)))
    bias.append(e/(2**n))

plt.plot(bias,sens,"o",markersize=0.2,label="ncf")
plt.xlabel("Activity bias")
plt.ylabel("Sensitivity")
plt.legend()
plt.savefig("partial_canalization.png")
plt.show()"""