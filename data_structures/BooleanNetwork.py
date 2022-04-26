from sample_functions import random_function
import numpy as np
import random
import time
import itertools


class BooleanNetwork:
    def __init__(self,node_number,in_degrees,random_method):
        if random_method==None:
            self.gene_functions=[]
        else:
            self.gene_functions=[random_function(in_degrees[i],random_method) for i in range(0,node_number)]
        connectivity_matrix=np.zeros((node_number,node_number))
        for i in range(0,node_number):
            row=[False]*(node_number-in_degrees[i]-1)+[True]*in_degrees[i]
            random.shuffle(row)
            row.insert(i,False)
            connectivity_matrix[i,:]=row
        self.connectivity_matrix=connectivity_matrix
        self.number_of_nodes=node_number
        self.in_degrees=in_degrees
        self.method=random_method
        self.sensitivity=None
        self.fixed_points={}
        self.cycles={}
    def update_gene(self,index,state):
        f=self.gene_functions[index]
        one_indices=np.where(self.connectivity_matrix[index,:]==1)[0]
        related_state=[state[i] for i in one_indices]
        return f[tuple(related_state)]*1

    def update_genes(self,state):
        return [self.update_gene(i,state) for i in range(0,self.number_of_nodes)]

    def find_sensitivity(self,index):
        logical_dict=self.gene_functions[index]
        if isinstance(logical_dict,dict):
            if not isinstance(logical_dict, dict) or logical_dict == {}:
                return 0, 0
            number_of_inputs = len(list(logical_dict.keys())[0])
            truth = [True, False]
            real_size = number_of_inputs
            s = 0
            for i in range(0, number_of_inputs):
                count = 0
                for c in itertools.product(truth, repeat=number_of_inputs - 1):
                    key = list(c)
                    key.insert(i, False)
                    temp = logical_dict[tuple(key)]
                    key[i] = True
                    if temp != logical_dict[tuple(key)]:
                        count += 1
                if count == 0:
                    real_size -= 1
                s += count / (2 ** (number_of_inputs - 1))
            return s, real_size
    def get_network_sensitivity(self):
        t=[]
        for i in range(self.number_of_nodes):
            t.append(self.find_sensitivity(i)[0])
        return sum(t)/len(t)
    def get_network_n_sensitivity(self,n):
        iterations=10000
        sum_hdist=0
        for i in range(iterations):
            state1 = [random.randint(0, 1) for i in range(self.number_of_nodes)]
            random.shuffle(state1)
            state2 = state1.copy()
            ri=random.randint(0,self.number_of_nodes-1)
            state2[ri]=1-state2[ri]
            for epoch in range(n):
                state1=self.update_genes(state1)
                state2=self.update_genes(state2)
            hdist=sum(np.array(state1) != np.array(state2))
            sum_hdist += hdist
        return sum_hdist/iterations

    def in_fixed_point(self,states):
        if states[-1]==states[-2]:
            try:
                self.fixed_points[tuple(states[-1])]+=1
            except:
                self.fixed_points[tuple(states[-1])]=1
            return True
        return False

    def in_cycle(self,states):
        if len(states)<2:
            return False
        temp=states.copy()
        last=temp.pop()
        try:
            i=temp.index(last)
        except:
            return False
        cycle = temp[i:]
        if len(cycle) == 1:
            return False
        t = is_new_cycle(cycle, self.cycles)
        if t == True:
            self.cycles[tuple(cycle)] = 1

        else:
            self.cycles[t] += 1
        return True

    def in_steady_state(self,states):
        if self.in_fixed_point(states):
            return True
        return self.in_cycle(states)

    def find_attractors(self,iteration=10000):
        self.cycles={}
        self.fixed_points={}
        for i in range(iteration):
            initial_state=[random.randint(0,1) for i in range(self.number_of_nodes)]
            random.shuffle(initial_state)
            states=[tuple(initial_state)]
            current_state=tuple(initial_state)
            while len(states)<2 or not self.in_steady_state(states):
                current_state=tuple(self.update_genes(current_state))
                states.append(current_state)

    def get_basin_entropy(self):
        if self.fixed_points or self.cycles:
            all_values=np.array(list(self.fixed_points.values())+list(self.cycles.values()))
            normalization=sum(all_values)
            normalized=all_values/normalization
            return sum(np.multiply(-normalized,np.log(normalized)))

def is_new_cycle(cycle, Cycles):
    for j in range(len(list(Cycles.keys()))):
        cand = list(Cycles.keys())[j]
        for i in range(0,len(cand)):
            temp = list(cand[1:])
            temp.append(cand[0])
            cand = temp
            if temp == cycle:
                return list(Cycles.keys())[j]
    return True

def gene_network_to_boolean(gene_reg):
    in_degrees=[gene.in_degree for gene in gene_reg.genes]
    bn=BooleanNetwork(len(gene_reg.genes),in_degrees,None)
    bn.gene_functions=[gene.construct_logical_dict() for gene in gene_reg.genes]
    return bn

