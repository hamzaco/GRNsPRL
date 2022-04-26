import itertools
import random
import matplotlib.pyplot as plt
import statistics
import json
import os
import math
import numpy as np
import pandas as pd
import altair as alt
import collections
from data_structures.Logic_Table_Network import logical_gene,random_symmetric_gene
from scipy.stats import pearsonr
plt.figure()
max_input=10
sample_size=1000
count2n=0
got_it=0
type0={'Sensitivity':[],'Bias':[],'input_size':[]}
type1={'Sensitivity':[],'Bias':[],'Distribution':[],'input_size':[],'small_class':[],'imbalance':[]}
type2={'Sensitivity':[],'Bias':[],'Distribution':[],'input_size':[],'small_class':[],'imbalance':[]}
type3={'Sensitivity':[],'Bias':[],'Distribution':[],'input_size':[],'small_class':[],'imbalance':[],'Sensitivity_appx':[]}

for n in np.arange(3,max_input+0.1,1):
    for a in np.arange(1,n+0.1,1):
        if a>n/2:
            continue
        minimum=min(a,n-a)
        maximum=max(a,n-a)
        imbalance=minimum/maximum
        s=-2 * n / 2 ** n + 2*a/2**a+ 2 * (n - a) / 2 ** (n - a)
        bias = (2 ** a - 1) * (2 ** (n - a) - 1) / (2 ** n)
        type3['Bias'].append(bias)
        type3['Bias'].append(1 - bias)
        k= 1/(2**a - 1)
        m= 1/ 2**a
        bias=1-bias
        s_appx= 2*a*m + 2*(-math.log2(k*(bias-m)))*(bias-m)-2*a*2**a*k*(bias-m)

        type3['Sensitivity_appx'].append(s_appx)
        type3['Sensitivity_appx'].append(s_appx)

        type3['Sensitivity'].append(s)
        type3['Sensitivity'].append(s)
        type3['small_class'].append(minimum)
        type3['small_class'].append(minimum)
        type3["imbalance"].append(imbalance)
        type3["imbalance"].append(imbalance)
        type3['input_size'].append(n)
        type3['input_size'].append(n)

        type3['Distribution'].append((a, n - a))
        type3['Distribution'].append((a, n - a))

        s = 2 * (n - 2 * a) / (2 ** n) + 2 * a / 2 ** a
        bias = (2 ** (n - a) - 1) / (2 ** n)
        type2['Sensitivity'].append(s)
        type2['Sensitivity'].append(s)
        type2['small_class'].append(minimum)
        type2['small_class'].append(minimum)
        type2["imbalance"].append(imbalance)
        type2["imbalance"].append(imbalance)

        type2['input_size'].append(n)
        type2['input_size'].append(n)

        type2['Distribution'].append((a, n - a))
        type2['Distribution'].append((a, n - a))
        type2['Bias'].append(bias)
        type2['Bias'].append(1 - bias)
        s = -2 * (n - 2 * a) / (2 ** n) + 2 * (n - a) / 2 ** (n - a)
        bias = (2 ** (a) - 1) / (2 ** n)
        type1['Sensitivity'].append(s)
        type1['Sensitivity'].append(s)
        type1['small_class'].append(minimum)
        type1['small_class'].append(minimum)
        type1["imbalance"].append(imbalance)
        type1["imbalance"].append(imbalance)

        type1['input_size'].append(n)
        type1['input_size'].append(n)

        type1['Distribution'].append((a, n - a))
        type1['Distribution'].append((a, n - a))
        type1['Bias'].append(bias)
        type1['Bias'].append(1 - bias)

        s= 2*n/2**n
        bias= 1/2**n
        type0['Sensitivity'].append(s)
        type0['Sensitivity'].append(s)
        type0['Bias'].append(bias)
        type0['Bias'].append(1 - bias)
        type0['input_size'].append(n)
        type0['input_size'].append(n)
"""plt.figure()
plt.plot(biases,sensitivities,'*k')
#plt.plot(biases_middle,sensitivities_middle,'*r')
plt.show()"""

print(statistics.mean(type1['Sensitivity']))
print(statistics.mean(type2['Sensitivity']))
print(statistics.mean(type3['Sensitivity']))
df0=pd.DataFrame(data=type0)
df1=pd.DataFrame(data=type1)
df2=pd.DataFrame(data=type2)
df3=pd.DataFrame(data=type3)

chart=alt.Chart(df3).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity',
    color='input_size',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','small_class','imbalance']
).interactive()

source = {'Sensitivity':[2*i/2**i for i in range(1,5)]*2, 'Bias':[1/2**i for i in range(1,5)]}
source['Bias'].extend([1-1/2**i for i in range(1,5)])
chart_points= alt.Chart(pd.DataFrame(source)).mark_point(size=60).encode(
    x='Bias',
    y='Sensitivity',
    color=alt.ColorValue('red')
)
alt.layer(chart,chart_points).show()
chart1=alt.Chart(df3).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity_appx',
    color='input_size',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','small_class','imbalance']
).interactive()


chart2=alt.Chart(df2).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity',
    color='input_size',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','small_class','imbalance']
).interactive()


alt.layer(chart2,chart_points).show()


chart3=alt.Chart(df1).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity',
    color='input_size',
    tooltip=['Sensitivity', 'Bias','Distribution','input_size','small_class','imbalance']
).interactive()

chart0=alt.Chart(df0).mark_circle(size=20).encode(
    x='Bias',
    y='Sensitivity',
    color='input_size',
    tooltip=['Sensitivity', 'Bias','input_size']
).interactive()

chart3.show()
chart0.show()

alt.layer(chart0,chart,chart2,chart3).show()