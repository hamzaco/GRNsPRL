
from data_structures.Logic_Table_Network import compare_logic_dict, restrict
import json
import os
import numpy as np
import statistics
import matplotlib.pyplot as plt
from data_structures.Gene_Network import Gene_Network, shuffle_network,print_genes
from data_structures.Gene import is_symmetric
import itertools
import math


def get_partitions(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in get_partitions(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset
        yield [ [ first ] ] + smaller
def get_layers(n,r):
    temp=set()
    for p in get_partitions([1]*n):
        p=[len(i) for i in p]
        if len(p)==r and p[-1]>=2:
            temp.add(tuple(p))
    return temp
def nCr(n,r):
    return math.factorial(n)/(math.factorial(n-r)*math.factorial(r))
def get_t_partition(s,k_list):
    temp=set()
    r=len(k_list)
    for p in get_partitions([1]*s):
        p=[len(i) for i in p]
        if len(p)==r and (3 not in p) and (sum(np.array(k_list) >= np.array(p))==r):
            temp.add(tuple(p))
    return temp

def multiplicity(n,k_list):
    denominator=1
    for k in k_list:
        denominator *= math.factorial(k)
    return math.factorial(n)/denominator
number_of_inputs=[3,4,5,6]
dim=max(number_of_inputs)+1
Nrs=np.zeros((dim,dim,dim))
Ns_ncf=np.zeros((dim,dim))
Nr=np.zeros((dim,dim))

for n in number_of_inputs:
    N=0
    for s in range(1,n):
        for r in range(math.ceil(s/2),s+1):
            sum_layers=0
            for layer_partition in get_layers(n,r):
                temp_sum = 0
                for t_partition in get_t_partition(s,layer_partition):
                    prod_t = 1
                    for i in range(len(t_partition)):
                        t=t_partition[i]
                        k=layer_partition[i]
                        prod_t *= ((t-1)*(2**k-2)+1-(-1)**t)
                    temp_sum += prod_t
                temp_sum *= multiplicity(n,layer_partition)
                sum_layers += temp_sum
            Nrs[n,r,s]= 2*sum_layers
            N += Nrs[n,r,s]
            Nr[n,r] += Nrs[n,r,s]
            print("n=",n,"\t s=",s,"\t r=",r,"\t N=", Nrs[n,r,s])
        Ns_ncf[n,s]=N
    temp=0
    for k in range(n-1):
        if n-2-k>=k:
            temp += nCr(n-2-k,k)*2**(n-2-2*k)
    Ns_ncf[n,n] = 2*math.factorial(n)*temp

