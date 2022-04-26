from data_structures.Gene_Network import Gene_Network, calculate_effect, search_gene
from data_structures import Gene
import itertools
import os
import json
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import pearsonr

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
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
with open('sensitivity_data/sens_vs_bias.json', 'r') as fp:
    sensitivity_vs_acitivity_bias = json.load(fp)
with open('dist_sens_bias_data/dist_sens_bias_brute_force.json', 'r') as fp:
    data_brute_force = json.load(fp)
with open('dist_sens_bias_data/dist_sens_bias_asymmetry.json', 'r') as fp:
    data_asymmetric = json.load(fp)
with open('dist_sens_bias_data/dist_sens_bias_symmetry_dictated.json', 'r') as fp:
    data_symmetry_dictated = json.load(fp)
with open('dist_sens_bias_data/dist_sens_bias_xor_excluded.json', 'r') as fp:
    data_xor_excluded = json.load(fp)
with open('dist_sens_bias_data/dist_sens_bias_canalizing.json', 'r') as fp:
    data_canalizing = json.load(fp)
#
#
#
#
#


#
#
#
#
#
sensitivity=[]
bias=[]
for input_number in range(2,len(data_canalizing)):
    for data in data_canalizing[input_number]:
        sensitivity.append(data[1])
        bias.append(data[2])

fig= plt.figure()
ax1 = fig.add_subplot(221,projection='3d')
hist, xedges, yedges = np.histogram2d(bias, sensitivity, bins=(50,50), range = [[0,1],[0,2]]) # you can change your bins, and the range on which to take data
xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:],indexing='ij') -(xedges[1]-xedges[0])

xpos = xpos.flatten()*1./2
ypos = ypos.flatten()*1./2
zpos = np.zeros_like(xpos)
dx = xedges [1] - xedges [0]
dy = yedges [1] - yedges [0]
dz = hist.ravel()
for i in range(len(dz)):
    if dz[i]>0:
        dz[i]=np.log(dz[i])
cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(dz)   # get range of colorbars so we can normalize
min_height = np.min(dz)
# scale each z to [0,1], and get their rgb values
rgba = [cmap((k-min_height)/max_height) for k in dz]
    

ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
plt.title("canalizing")


#
#
#
#
#
#
#
#
#
#
sensitivity = []
bias = []
for input_number in range(2, len(data_xor_excluded)):
    for data in data_xor_excluded[input_number]:
        sensitivity.append(data[1])
        bias.append(data[2])

ax2 = fig.add_subplot(222,projection='3d')
hist, xedges, yedges = np.histogram2d(bias, sensitivity, bins=(50, 50), range=[[0, 1], [0,
                                                                                          2]])  # you can change your bins, and the range on which to take data
xpos, ypos = np.meshgrid(xedges[:-1] + xedges[1:], yedges[:-1] + yedges[1:], indexing='ij') - (xedges[1] - xedges[0])

xpos = xpos.flatten() * 1. / 2
ypos = ypos.flatten() * 1. / 2
zpos = np.zeros_like(xpos)
dx = xedges[1] - xedges[0]
dy = yedges[1] - yedges[0]
dz = hist.ravel()
for i in range(len(dz)):
    if dz[i] > 0:
        dz[i] = np.log(dz[i])
cmap = plt.cm.get_cmap('jet')  # Get desired colormap - you can change this!
max_height = np.max(dz)  # get range of colorbars so we can normalize
min_height = np.min(dz)
# scale each z to [0,1], and get their rgb values
rgba = [cmap((k - min_height) / max_height) for k in dz]

ax2.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')
plt.title("xor excluded")


#
#
#
#
#
#
#
#
#
#
sensitivity=[]
bias=[]
for data in sensitivity_vs_acitivity_bias:
    sensitivity.append(data[0])
    bias.append(data[1])

ax3 = fig.add_subplot(223,projection='3d')
hist, xedges, yedges = np.histogram2d(bias, sensitivity, bins=(50,50), range = [[0,1],[0,2]]) # you can change your bins, and the range on which to take data
xpos, ypos = np.meshgrid(xedges[:-1]+xedges[1:], yedges[:-1]+yedges[1:],indexing='ij') -(xedges[1]-xedges[0])

xpos = xpos.flatten()*1./2
ypos = ypos.flatten()*1./2
zpos = np.zeros_like(xpos)
dx = xedges [1] - xedges [0]
dy = yedges [1] - yedges [0]
dz = hist.ravel()
for i in range(len(dz)):
    if dz[i]>0:
        dz[i]=np.log(dz[i])
cmap = plt.cm.get_cmap('jet') # Get desired colormap - you can change this!
max_height = np.max(dz)   # get range of colorbars so we can normalize
min_height = np.min(dz)
# scale each z to [0,1], and get their rgb values
rgba = [cmap((k-min_height)/max_height) for k in dz]

ax3.bar3d(xpos, ypos, zpos, dx, dy, dz, color=rgba, zsort='average')

#ax3.view_init(90, 0)
plt.title("Original")
plt.savefig('2D_histrograms_sens_vs_bias_50bins')
plt.show()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
"""sensitivity=[]
bias=[]
for input_number in range(2,len(data_symmetry_dictated)):
    for data in data_symmetry_dictated[input_number]:
        sensitivity.append(data[1])
        bias.append(data[2])
plt.figure()
plt.hist2d(bias,sensitivity,bins=[1,2])
plt.show()"""