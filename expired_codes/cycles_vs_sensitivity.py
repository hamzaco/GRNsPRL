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

with open('sensitivities_attractors.json', 'r') as fp:
    sensitivities = json.load(fp)

with open('attractors.json', 'r') as fp:
    attractors = json.load(fp)


cycle_percentage={}
plt.figure()
i=0
cps=np.zeros(len(sensitivities))
s=np.zeros(len(sensitivities))
for cell in sensitivities:
    att=attractors[cell]
    cycles=0
    for attractor in att:
        if not isinstance(attractor[0],str):
            cycles +=1
    cycle_percentage[cell]=cycles/len(att)

    plt.plot(cycle_percentage[cell],sensitivities[cell],'*r')
    print(sensitivities[cell])
    if cycles==0:
        continue
    cps[i]=cycle_percentage[cell]
    s[i]=sensitivities[cell]
    i+=1
plt.plot(cycle_percentage['CC_Arabidopsis_thaliana_CC'],sensitivities['CC_Arabidopsis_thaliana_CC'],'*r',label='Original')
for i in range(0,len(cps)):
    if cps[i]==0:
        s = np.delete(s, np.arange(i, len(cps)))
        cps=np.delete(cps,np.arange(i,len(cps)))

        break
p=np.polyfit(cps,s,1)
appx_s=np.polyval(p,cps)

plt.plot(cps,appx_s,label='Fit')
plt.legend()
plt.ylabel('Sensitivities')
plt.xlabel('Cycle Percentage')
plt.savefig('Cycle Percentage vs Sensitivity')
plt.show()
