import matplotlib.pyplot as plt
import numpy as np
import json


main_dir="C:/Users/Lenovo/PycharmProjects/INDEP"
with open(main_dir+'/general-data/average_inputs.json', 'r') as fp:
    average_inputs = json.load(fp)
with open(main_dir+'/sensitivity_data/sensitivities.json', 'r') as fp:
    sensitivities = json.load(fp)
with open(main_dir+'/sensitivity_data/sensitivities_random_symmetric.json', 'r') as fp:
    sensitivities_random = json.load(fp)
with open(main_dir+'/error_bars/error_bars_random_sensitivities_symmetric.json', 'r') as fp:
    error_bars = json.load(fp)
with open(main_dir+'/sensitivity_data/sensitivities_totally_random_.json', 'r') as fp:
    sensitivities_totally_random = json.load(fp)

temp=list(sensitivities.values())
temp = [tup[0] for tup in temp]
plt.figure()
plt.hist(temp, bins=30, range=[0,3], facecolor='blue', ec='black', alpha=0.2, label='Original')

entries, edges, _ = plt.hist(list(sensitivities_random.values()), bins=30, range=[0,3], facecolor='red', ec='black', alpha=0.2, label='Random (Symmetric)')
bin_centers = 0.5 * (edges[:-1] + edges[1:])
plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='r_',alpha=0.2)

entries, edges, _ = plt.hist(list(sensitivities_totally_random.values()), bins=30, range=[0,3], facecolor='green', ec='black', alpha=0.2, label='Totally Random')
bin_centers = 0.5 * (edges[:-1] + edges[1:])
plt.errorbar(bin_centers, entries, yerr=np.sqrt(entries), fmt='g_',alpha=0.2)

plt.legend()
plt.xlabel('Sensitivity')
plt.ylabel('Count')

#plt.savefig('Histogram sensitivities')
plt.show()