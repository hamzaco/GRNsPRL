from data_structures.Gene_Network import Gene_Network, calculate_effect, search_gene
import json
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
from pylab import cm

main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"

path = "C:/Users/Lenovo/PycharmProjects/INDEP/network_data_altered"
path_cells = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"

with open(main_dir+'/sensitivity_data/sens_vs_bias.json', 'r') as fp:
    sensitivity_vs_acitivity_bias = json.load(fp)

with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_ncf.json', 'r') as fp:
    data_ncf = json.load(fp)

with open(main_dir + '/sens_bias_stdev_total_random_cell_collective.txt', 'r') as fp:
    sens_bias_stdev_random = json.load(fp)

sens_bias_stdev_random=sens_bias_stdev_random
low_line=[data[0]-data[1] for data in sens_bias_stdev_random]
middle_line=[data[0] for data in sens_bias_stdev_random]
high_line=[data[0]+data[1] for data in sens_bias_stdev_random]
bias_random=np.arange(0,1.0,step=0.07)
print(bias_random)
print(high_line)

print(low_line)

biases=[]
sensitivities=[]
choosen_data= data_ncf
for input_number in range(2,len(choosen_data)):
    for data in choosen_data[input_number]:
        bias = data[2]
        sensitivities.append(data[1])
        biases.append(bias)

ind=np.argsort(biases)

biases_ncf=np.take_along_axis(np.array(biases),ind,axis=0)
sensitivities_ncf=np.take_along_axis(np.array(sensitivities),ind,axis=0)

biases=[]
sensitivities=[]
biases_hist=[]
sensitivities_hist=[]
freq_table={}

for data in sensitivity_vs_acitivity_bias:
    biases_hist.append(data[1])
    sensitivities_hist.append(data[0])
    try:
        freq_table[tuple(data)]+=1
    except:
        freq_table[tuple(data)]=1

print(len(sensitivities_hist))
freq_table= dict(sorted(freq_table.items(), key=lambda item: item[1]))

for key in freq_table.keys():
    sensitivities.append(key[0])
    biases.append(key[1])
colors=list(map(lambda x: np.log(x),list(freq_table.values())))



def scatter_hist(x1, y1,x2,y2, x3,y3,ax, ax_histx, ax_histy):
    # no labels
    ax_histx.tick_params(axis="x",labelsize=15)
    ax_histx.tick_params(labelsize=15)
    ax_histy.tick_params(labelsize=15)
    ax_histy.tick_params(axis="y",labelsize=15)
    ax.tick_params(labelsize=15)

    # the scatter plot:


    ax.plot(x2, y2,label="NCFs",color="k")


    colors = list(map(lambda x: np.log(x), list(freq_table.values())))
    ax.scatter(x1, y1,s=20,c=colors, cmap="rainbow",marker="o",label="CC Database")
    cbar=fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,vmax=max(list(freq_table.values()))/sum(list(freq_table.values()))), cmap="rainbow"))
    cbar.ax.tick_params(labelsize=14)

    # now determine nice limits by hand:
    binwidth = 0.05
    ax.fill_between(bias_random,low_line,high_line,alpha=0.05,color="k",label="Random Boolean functions")
    ax.plot([-0.5,1.5],[1,1],c="red",alpha=0.5,label="Edge of chaos")
    aratio=np.arange(0,1.01,step=0.05)
    #ax.plot(aratio,4*aratio*(1-aratio),"--",alpha=0.5,label="Graph theory lower bound")
    ax.set_ylabel("Sensitivity "+r'($\xi$)',fontsize=20)
    ax.set_xlabel("Activity Ratio "+r'($p$)',fontsize=20)

    bins_sens = np.arange(binwidth/2, 2, step=binwidth)
    bins_bias= np.arange(-binwidth/4,1.05,step=binwidth/2)
    ax_histx.hist(x3, bins=bins_bias,ec="k")
    ax_histy.hist(y3, bins=bins_sens,ec="k", orientation='horizontal')

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005


rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

# start with a square Figure
fig = plt.figure(figsize=(11, 10))

ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function

scatter_hist(biases, sensitivities, biases_ncf,sensitivities_ncf,
             biases_hist,sensitivities_hist,ax, ax_histx, ax_histy)
ax.set_xlim(-0.05,1.05)
ax.set_ylim(-0.05,2.55)
ax.set_xticks(np.arange(0,1.1,step=0.1))
ax.set_yticks(np.arange(0,2.51,step=0.2))
#get handles and labels
handles, labels = ax.get_legend_handles_labels()
# sort both labels and handles by labels
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
print(labels)

order=[0,1,2,3]
labels=[labels[idx] for idx in order]
handles=[handles[idx] for idx in order]
ax.legend(handles, labels,loc=(0.24,0.02),fontsize=15)
#ax.legend()
#plt.figtext(0.865,0.78,"Log scale",fontsize=13)
plt.figtext(0.92,0.76,"Freq.",fontsize=18)
plt.figtext(0.11,0.37,"Unstable",fontsize=16)
plt.figtext(0.11,0.34,"Stable",fontsize=16)

plt.savefig("sens_vs_bias_original_with_histogram_2_included_equal_prob_simple.png")
plt.show()

plt.plot(biases_ncf,sensitivities_ncf,c="k")
plt.xticks(np.arange(0,1.1,step=0.1),fontsize=11)
plt.yticks(np.arange(0,1.5,step=0.2),fontsize=11)
plt.xlabel("Activity ratio",fontsize=16)
plt.ylabel("Sensitivity",fontsize=16)
plt.savefig("blanchmange_sens_bias.png")
plt.show()

plt.plot(biases_ncf,sensitivities_ncf,c="k")
plt.axis("off")
plt.savefig("blanchmange_sens_bias_no_axis.png")
plt.show()