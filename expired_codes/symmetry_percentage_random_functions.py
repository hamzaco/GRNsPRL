import itertools
import random
import matplotlib.pyplot as plt
import statistics
import json
import os
from data_structures.Gene import Gene
import math
from Logic_Table_Network import logical_gene
from scipy.stats import pearsonr
truth=[True,False]
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)
max_input=8
iteration_size=5000
with open('inputs_and_sensitivities.json','r') as fp:
    inputs_and_sensitivities=json.load(fp)

activity_bias=[]
no_symmetry_activity_bias=[]
truth=[True,False]
no_symmetry_count_vs_input=[]
no_symmetry_sensitivities=[]
sensitivities=[]
plt.figure()
for input in range(2,max_input):
    odd_numbers=[]
    for i in range(1,2**(input)):
        odd_numbers.append(i)
    count = 0
    all=0
    for j in range (0,iteration_size):
        logical_dict={}
        number_of_truths=random.choice(odd_numbers)
        truths = [True] * number_of_truths + [False] * (2 **input - number_of_truths)

        random.shuffle(truths)
        i=0
        for c in itertools.product(truth, repeat=input):
            logical_dict[c]=random.choice(truths)
            i+=1
        gene=Gene(name='random '+str(j), function=logical_dict)

        s,real_input_size=gene.find_sensitivity()
        bias=number_of_truths/(2**input)

        if real_input_size != input:
            j = j - 1
            continue
        all +=1
        dist, classes = gene.find_symmetries()
        if dist== [1]*input:
            count +=1
            no_symmetry_sensitivities.append(s)
            no_symmetry_activity_bias.append(bias)
        else:
            sensitivities.append(s)
            activity_bias.append(bias)
            plt.plot(bias,s,'*k')
    print(input,1-count/all)
    no_symmetry_count_vs_input.append(count)
plt.show()


plt.figure()

theoretical = [iteration_size*(1-2**(-2**(i-2)+1))**(nCr(i,2)) for i in range(2,max_input)]
plt.plot(no_symmetry_count_vs_input,'*k')
plt.plot(theoretical,'*r')
plt.show()

"""
plt.figure
plt.hist(sensitivities,bins=40,range=[0,4],weights=[1/len(sensitivities)]*len(sensitivities), facecolor='red', ec='black', alpha=0.2, label='Symmetry')
plt.hist(no_symmetry_sensitivities,bins=40,range=[0,4],weights=[1/len(no_symmetry_sensitivities)]*len(no_symmetry_sensitivities), facecolor='blue', ec='black', alpha=0.2, label='No symmetry')
plt.legend()
plt.xlabel('sensitivity')
#plt.savefig('random_nodes_sensitivies_symmetry_1-5 inputs')
plt.show()

plt.figure()
n_symmetric,bins, _= plt.hist(activity_bias, bins=20, range=[0, 1], facecolor='red',
         ec='black', alpha=0.2,weights=[1/len(activity_bias)]*len(activity_bias), label='Symmetry')
n_nosymmetry,_,_=plt.hist(no_symmetry_activity_bias, bins=20, range=[0, 1],facecolor='blue', ec='black',
         alpha=0.2,weights=[1/len(no_symmetry_activity_bias)]*len(no_symmetry_activity_bias),label='No symmetry')
n_original,_,_=plt.hist(activity_bias_original, bins=20, range=[0, 1], facecolor='green', ec='black',
         alpha=0.2,weights=[1/len(activity_bias_original)]*len(activity_bias_original) ,label='original')
plt.xlabel('Activity bias')
plt.legend()
#plt.savefig('random_nodes_activity_bias_symmetry_1-5 inputs')
plt.show()

plt.figure()


plt.loglog(n_original/sum(n_original),n_symmetric/(sum(n_symmetric)),'*r',label='Symmetric')
plt.loglog(n_original/sum(n_original),n_nosymmetry/(sum(n_nosymmetry)),'*g', label= 'No symmetry')
plt.loglog(n_original/sum(n_original),n_original/sum(n_original),'*k')
plt.grid(True)
plt.legend()
plt.ylabel('Random symmetric activity bias')
plt.xlabel('Original activity bias')
#plt.savefig('Acitivity bias(randoms vs original)')
plt.show()

corr, _ = pearsonr(n_symmetric,n_original)
print(corr)
corr, _ = pearsonr(n_nosymmetry,n_original)
print(corr)"""