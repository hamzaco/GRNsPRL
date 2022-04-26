import pandas as pd
import altair as alt
import json
from statistics import mean

def get_type(data):
    a=min(data[0])
    n=sum(data[0])
    s=data[1]
    bias=data[2]
    bias= 1/2- abs(bias-1/2)
    if bias == (2**a-1)/2**n and s==2*(2*a-n)/2**n + 2*(n-a)/2**(n-a):
        return 1
    if bias == 1/2**a-1/2**n and s==2*a/2**a + (2*n-4*a)/2**n :
        return 2
    if bias ==  1/2**a+(2**a-1)/2**n and s== -2 * n / 2 ** n + 2*a/2**a+ 2 * (n - a) / 2 ** (n - a):
        return 3
    return 4
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original.json', 'r') as fp:
    data_original = json.load(fp)
with open(main_dir+'/dist_sens_bias_data/dist_sens_bias_original_3_or_more.json', 'r') as fp:
    data_original_large = json.load(fp)
two_class_count=0
all_count=0
one_class_count=0
asymmetry_count=0
sensitivities_all=[]
type_data=[[]]*5
for input_number in range(3,len(data_original)):
    for data in data_original[input_number]:
        if data[1]>0:
            all_count += 1
            if len(data[0])==input_number:
                asymmetry_count += 1
            if len(data[0])==1:
                one_class_count +=1
                temp = type_data[0].copy()
                temp.append(data[1])
                type_data[0] = temp
            elif len(data[0])==2:
                temp=type_data[get_type(data)].copy()
                temp.append(data[1])
                type_data[get_type(data)]=temp
                two_class_count +=1
                sensitivities_all.append(data[1])
dist_type= [len(i) for i in type_data]
print(dist_type)
mean_type= [mean(i) for i in type_data]
all_mean = sum([sum(i) for i in type_data])/(one_class_count+two_class_count)
print(mean_type)
print(all_mean)
print("All genes: "+str(all_count))
print("Asymmetric genes: "+str(asymmetry_count))
print("One class percentage: "+str(one_class_count/(all_count)))
print("At most two class percentage: "+str((one_class_count+two_class_count)/all_count))
print("Two class percentage: "+str(two_class_count/(all_count-one_class_count)))
sensitivities_large=[]
for input_number in range(0,len(data_original_large)):
    for data in data_original_large[input_number]:
        if data[1]>0:
            if len(data[0])>2:
                sensitivities_large.append(data[1])

print(mean(sensitivities_large))
alt.data_transformers.disable_max_rows()
sensitivities={'Sensitivities_all':sensitivities_all}
sensitivities=pd.DataFrame(sensitivities)
chart1=alt.Chart(sensitivities).transform_joinaggregate(
    total='count(*)'
).transform_calculate(
    pct='1 / datum.total'
).mark_bar(
).encode(
    alt.X("Sensitivities_all:Q",  bin=alt.Bin(step=0.01)),
    alt.Y('sum(pct):Q', axis=alt.Axis(format='%')),
)

sensitivities={'Sensitivities_large':sensitivities_large}
sensitivities=pd.DataFrame(sensitivities)
chart2=alt.Chart(sensitivities).transform_fold(
    ['Sensitivities_large '],
).transform_joinaggregate(
    total='count(*)'
).transform_calculate(
    pct='1 / datum.total'
).mark_bar(
).encode(
    alt.X("Sensitivities_large:Q", bin=alt.Bin(step=0.1)),
    alt.Y('sum(pct):Q', axis=alt.Axis(format='%')),
)

all=alt.vconcat(chart1,chart2)

all.show()
