import itertools
import random
import matplotlib.pyplot as plt
import statistics
import json
import os
import math
import pandas as pd
from data_structures.Logic_Table_Network import logical_gene
from scipy.stats import pearsonr
import altair as alt
type1={'Sensitivity':[],'Bias':[],'Distribution':[],'input_size':[],'imbalance':[],"small_class":[],"large_class":[]}

near_one=0
all=0
near_one_2_classes=0
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original.json','r') as fp:
    dist_sens_bias=json.load(fp)
imbalances=[]
near_one_small_classes=[]
near_one_len_dist=[]
far_one_small_classes=[]
for input_number in range(2,len(dist_sens_bias)):
    datas=dist_sens_bias[input_number]
    for data in datas:
        dist = data[0]
        all +=1
        if data[1]<=1.05 and data[1]>=0.95:
            near_one +=1
            near_one_len_dist.append(len(dist))
            near_one_small_classes.append(min(dist))
        elif data[1]>1.20 or data[1]<0.80:
            far_one_small_classes.append(min(dist))
        else:
            continue
        type1["Sensitivity"].append(data[1])
        type1["small_class"].append(min(dist))
        type1["large_class"].append(max(dist))

        type1["Bias"].append(data[2])
        type1["Distribution"].append(dist)
        type1["input_size"].append(input_number)
        type1["imbalance"].append(1)

df=pd.DataFrame(type1)

plt.hist(near_one_len_dist,range=[1,5])
plt.show()
plt.hist(near_one_small_classes,range=[1,6])
plt.show()
plt.hist(far_one_small_classes,range=[1,6])
plt.show()

chart=alt.Chart(df).mark_circle(size=20).encode(
    x='small_class',
    y='large_class',
    color='imbalance',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','imbalance']
).interactive()

chart.show()