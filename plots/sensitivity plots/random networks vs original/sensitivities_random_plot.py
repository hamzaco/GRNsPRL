import matplotlib.pyplot as plt
import numpy as np
import json
main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"

with open(main_dir+'/general-data/average_inputs.json', 'r') as fp:
    average_inputs = json.load(fp)
with open(main_dir+'/sensitivity_data/sensitivities.json', 'r') as fp:
    sensitivities = json.load(fp)
with open(main_dir+'/sensitivity_data/sensitivities_random_brute_force.json', 'r') as fp:
    sensitivities_random = json.load(fp)
with open(main_dir+'/error_bars/error_bars_random_odd_sensitivities_symmetric.json', 'r') as fp:
    error_bars = json.load(fp)
with open(main_dir+'/p_values/p_values_random_odd_sensitivities_symmetric.json', 'r') as fp:
    p_values = json.load(fp)
big_entries= [
"CC_Inflammatory_Bowel_Disease_IBD_Model",
"CC_Yeast_Apoptosis"]

y= np.zeros(len(sensitivities_random))
yrand =  np.zeros(len(sensitivities_random))
errors= np.zeros(len(sensitivities_random))
x= np.arange(0,len(sensitivities_random))
i=0
t=[]
for cell in error_bars:
    if cell in big_entries:
        continue

    y[i]=sensitivities[cell][0]
    yrand[i] = sensitivities_random[cell]
    errors[i] = error_bars[cell]
    t.append((yrand[i],y[i],errors[i],average_inputs[cell]))
    i +=1
t.sort()
y.sort()
yrand.sort()
plt.figure()
plt.plot(x[0],y[0],'*r',label='original')
plt.plot(x[0], yrand[0], '*k', label='random' )
for i in range(0,len(y)):
    temp=t[i]
    plt.plot(x[i],yrand[i],'*k')
    plt.plot(x[i], y[i], '*r' )
    plt.errorbar(x[i], yrand[i], xerr=None, yerr=temp[2], fmt='b')
plt.xlabel('Network number')
plt.ylabel('Sensitivities')
plt.grid(True)
plt.legend()
plt.title('Sensitivities of Original and Random (Odd and Symmetric) Networks ')
#plt.savefig('sensitivities_random_odd_symmetric')
plt.show()

