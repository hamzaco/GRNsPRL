from Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
path = "C:/Users/Lenovo/PycharmProjects/INDEP"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

big_entries= ["CC_CD4p_Tcell_Diff",
"CC_Differentiation_T_lymphocytes",
"CC_Inflammatory_Bowel_Disease_IBD_Model",
"CC_Glucose_Repression_Signaling_2009",
"CC_HCC1954_Breast_Cell_Line_Long-term_ErbB_Network",
"CC_IGVH_mutations",
"CC_Lymphopoiesis_Regulatory_Network",
"CC_Signaling_in_Macrophage_Activation",
"CC_Signal_Transduction_in_Fibroblasts",
"CC_Treatment_of_Castration_Resistant_Prostate_Cancer",
"CC_Yeast_Apoptosis"]
entries = os.listdir(path_cells)
with open(main_dir+'/sensitivity_data/sens_vs_bias.json', 'r') as fp:
    sensitivity_vs_acitivity_bias = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original.json', 'r') as fp:
    data_original = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_brute_force.json', 'r') as fp:
    data_brute_force = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_asymmetry.json', 'r') as fp:
    data_asymmetric = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_symmetry_dictated.json', 'r') as fp:
    data_symmetry_dictated = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_xor_excluded.json', 'r') as fp:
    data_xor_excluded = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_canalizing.json', 'r') as fp:
    data_canalizing = json.load(fp)

activity_bias_original=[]
brute_force_activity_bias=[]
symmetry_dictation_activity_bias=[]
asymmetry_activity_bias=[]
xor_excluded_activity_bias=[]
canalizing_activity_bias=[]
for input_number in range(2,len(data_brute_force)):
    for data in data_brute_force[input_number]:
        if data[2]>0.5:
            brute_force_activity_bias.append(1-data[2])
        else:
            brute_force_activity_bias.append(data[2])


for input_number in range(2,len(data_xor_excluded)):
    for data in data_symmetry_dictated[input_number]:
        if data[2] > 0.5:
            xor_excluded_activity_bias.append(1 - data[2])
        else:
            xor_excluded_activity_bias.append(data[2])

for input_number in range(2,len(data_asymmetric)):
    for data in data_asymmetric[input_number]:
        if data[1] > 0.5:
            asymmetry_activity_bias.append(1 - data[1])
        else:
            asymmetry_activity_bias.append(data[1])

for input_number in range(2,len(data_original)):
    for data in data_original[input_number]:
        if data[2]>0.5:
            activity_bias_original.append(1-data[2])
        else:
            activity_bias_original.append(data[2])

plt.figure()
n_symmetric,bins, _= plt.hist(xor_excluded_activity_bias, bins=10, range=[0, 0.5], facecolor='green',
         ec='black', alpha=0.2,weights=[1/len(xor_excluded_activity_bias)]*len(xor_excluded_activity_bias), label='Xor excluded')
n_brute_force,_,_=plt.hist(brute_force_activity_bias, bins=10, range=[0, 0.5],facecolor='blue', ec='black',
         alpha=0.2,weights=[1/len(brute_force_activity_bias)]*len(brute_force_activity_bias),label='Brute force')

n_original,_,_=plt.hist(activity_bias_original, bins=10, range=[0, 0.5], facecolor='red', ec='black',
         alpha=0.2,weights=[1/len(activity_bias_original)]*len(activity_bias_original) ,label='original')
n_asymmetry,_,_=plt.hist(asymmetry_activity_bias, bins=10, range=[0, 0.5], facecolor='orange', ec='black',
         alpha=0.2,weights=[1/len(asymmetry_activity_bias)]*len(asymmetry_activity_bias) ,label='asymmetry')

plt.xlabel('Activity bias')
plt.legend()
plt.show()

plt.figure()
plt.hist([brute_force_activity_bias,
          activity_bias_original],bins=10, range=[0, 0.5],
         weights=[[1/len(brute_force_activity_bias)]*len(brute_force_activity_bias),
                [1/len(activity_bias_original)]*len(activity_bias_original) ],
         label=['Brute Force','Original'])
plt.xlabel('Activity bias')
plt.legend()
plt.savefig('brute_force_vs_original_activity_bias')
plt.show()


plt.figure()
plt.hist([asymmetry_activity_bias,
          activity_bias_original],bins=10, range=[0, 0.5],
         weights=[[1/len(asymmetry_activity_bias)]*len(asymmetry_activity_bias),
                [1/len(activity_bias_original)]*len(activity_bias_original) ],
         label=['Asymmetric','Original'])
plt.xlabel('Activity bias')
plt.legend()
plt.savefig('asymmetry_vs_original_activity_bias')
plt.show()


plt.figure()
plt.hist([xor_excluded_activity_bias,
          activity_bias_original],bins=10, range=[0, 0.5],
         weights=[[1/len(xor_excluded_activity_bias)]*len(xor_excluded_activity_bias),
                [1/len(activity_bias_original)]*len(activity_bias_original) ],
         label=['Xor excluded','Original'])
plt.xlabel('Activity bias')
plt.legend()
plt.savefig('symmetry_dictation_vs_original_activity_bias')
plt.show()

plt.figure()


plt.loglog(n_original/sum(n_original),n_symmetric/(sum(n_symmetric)),'*r',label='Xor excluded')
plt.loglog(n_original/sum(n_original),n_brute_force/(sum(n_brute_force)),'*g', label= 'Brute force')
plt.loglog(n_original/sum(n_original),n_original/sum(n_original),'*k')
plt.loglog(n_original/sum(n_original),n_asymmetry/sum(n_asymmetry),'*m', label='Asymmetric')

plt.grid(True)
plt.legend()
plt.ylabel('Random symmetric activity bias')
plt.xlabel('Original activity bias')
plt.savefig('Acitivity bias(randoms vs original)')
plt.show()

n_symmetric= [np.log(i) for i in n_symmetric]
n_original= [np.log(i) for i in n_original]
n_brute_force= [np.log(i) for i in n_brute_force]
n_asymmetry= [np.log(i) for i in n_asymmetry]

corr, _ = pearsonr(n_symmetric,n_original)
print("Symetric correlation",corr)
corr, _ = pearsonr(n_brute_force,n_original)
print("Brute force correlation",corr)
corr, _ = pearsonr(n_asymmetry,n_original)
print("Asymetric correlation",corr)