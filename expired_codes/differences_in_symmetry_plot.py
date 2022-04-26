import matplotlib.pyplot as plt
import numpy as np
import json
import itertools

with open('differences_in_symmetry_totally random.json', 'r') as fp:
    differences_in_symmetry = json.load(fp)

cumulative_differences=[]
for cell in differences_in_symmetry:
    cum_diff=sum(differences_in_symmetry[cell])/len(differences_in_symmetry[cell])
    cumulative_differences.append(cum_diff)

plt.figure()
plt.hist(cumulative_differences,bins=20, range=[-1,1],weights=[1/len(differences_in_symmetry)]*len(cumulative_differences))
plt.xlabel('Percentage of differences in symmetry')
plt.ylabel('Probability')
plt.savefig('percentage_differences_in_symmetry_random')
plt.show()