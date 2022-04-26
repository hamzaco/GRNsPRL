
from Logic_Table_Network import logical_gene,Logical_Table_Network,random_network
import json
import os
import numpy as np
import scipy as sp
from simulation import simulate
import statistics
import matplotlib.pyplot as plt
from Gene_Network import Gene_Network, shuffle_network,print_genes
from scipy.stats import pearsonr
import itertools
with open('and_vs_sensitivity.json','r') as fp:
    and_vs_sens=json.load(fp)
for desired_input in range (3,15):
    plt.figure()
    and_percentages=[]
    for data in and_vs_sens:
        s=data[0]
        input_number=s[1]
        s=s[0]

        if input_number==desired_input:
            plt.plot(data[1],s,'*k')
            and_percentages.append(data[1])
    plt.title('And percentage vs. Sensitivity (' + str(desired_input)+ '-input)')
    plt.xlabel('And Percentage')
    plt.ylabel('Sensitivity')
    #plt.savefig('and_vs_sensitivity_'+str(desired_input)+'_input')
    #plt.show()
    if len(and_percentages) >0:
        plt.figure()
        plt.hist(and_percentages,bins=10,range=[0, 1],
                 weights=[1/len(and_percentages)]*len(and_percentages),facecolor='orange', ec='black')
        plt.xlabel('And percentages')
        plt.title('Distrubition of and percentages in original networks (' + str(desired_input)+ '-input)')
        plt.ylabel('Probability')
        plt.savefig('and_distribution_original _'+str(desired_input)+'_input')
        plt.show()