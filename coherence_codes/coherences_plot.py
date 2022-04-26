import matplotlib.pyplot as plt
import numpy as np
import json
with open('sensitivities_neet.json', 'r') as fp:
    coherences = json.load(fp)
with open('sensitivities_totally_random_.json', 'r') as fp:
    coherences_random = json.load(fp)
with open('error_bars_totally_random_sensitivities.json', 'r') as fp:
    error_bars = json.load(fp)
with open('p_values_totally_random_sensitivities.json', 'r') as fp:
    p_values = json.load(fp)
y= np.zeros(len(coherences_random))
yrand =  np.zeros(len(coherences_random))
errors= np.zeros(len(coherences_random))
x= np.arange(0,len(coherences_random))
i=0
t=[]
print(coherences_random)
for cell in coherences_random:
    y[i]=coherences[cell]
    yrand[i] = coherences_random[cell]
    errors[i] = error_bars[cell]
    t.append((y[i],yrand[i],errors[i]))
    i +=1
t.sort()
plt.figure()
for i in range(0,len(t)):
    temp=t[i]
    print(i),
    print("th error:"),
    print(temp[2])

    plt.plot(x[i],temp[0],'*r',label='original')
    plt.plot(x[i], temp[1], '*k', label='random' )
    plt.errorbar(x[i],temp[1],xerr=None,yerr=temp[2],fmt='b')
#plt.xlabel('Network number')
plt.ylabel('Sensitivity')
plt.title('Sensitivities Original vs. Totally Random')
plt.savefig('sensitivitiy_totally_random')
plt.show()

temp=list(p_values.values())
temp.sort()
plt.hist(temp,bins=10,range=[0,11])
plt.title('P-values of Original Coherence Coefficient with Random Distribution')
plt.savefig('p-values_sensitivity_totally_random')
plt.show()