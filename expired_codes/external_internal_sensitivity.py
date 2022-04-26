from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from Gene_Network import Gene_Network, shuffle_network,print_genes
import itertools

with open('sensitivities_external_seperated.json', 'r') as fp:
    sensitivities = json.load(fp)

datas=list(sensitivities.values())
datas.sort()

plt.figure()
i=1
c=0
for data in datas:
    if data[1]==0:
        continue
    plt.plot(i,data[0],'*r')
    plt.plot(i,data[1],'*k')
    if data[0]<data[1]:
        c+=1
    i+=1
plt.plot(i, data[0], '*r',label='Internal sensitivity')
plt.plot(i, data[1], '*k',label='External sensitivity')
plt.legend()
plt.ylabel('Sensitivities')
plt.show()
print(c)
print(i)