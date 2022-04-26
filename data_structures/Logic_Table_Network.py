from data_structures.Gene import Gene
import random
import itertools
import numpy as np
from data_structures.Gene_Network import Gene_Network,optimize_texts
class logical_gene:
    def __init__(self, name="", activation=False, function={}):
        self.name = name
        self.activation = activation #boolean
        self.logical_dict = function
        self.is_influencing = set() #set of genes
        self.is_influenced = set()  #set of genes
        self.is_influenced_externals = set() #set of strings
        self.negative_signals = set()
        self.positive_signals = set()
        self.number_of_inputs= None
        self.number_of_truths=0
    def activate(self):
        self.activation = True

    def deactivate(self):
        self.activation = False

    def find_structure(self):
        distribution, classes=self.find_symmetries()
        inputs = [i.name for i in list(self.is_influenced)]
        inputs = inputs + list(self.is_influenced_externals)
        positives = []
        for inpu in self.positive_signals:
            if isinstance(inpu,str):
                positives.append(inpu)
            else:
                positives.append(inpu.name)
        operation_relation=[]

        for pair in itertools.combinations(classes,2):
            Conditions = {}
            others=classes.copy()
            others.remove(pair[0])
            others.remove(pair[1])
            first=pair[0]
            second=pair[1]
            others_input=[]
            for inpu in inputs:
                others_input.append(inpu)

            for other_state in itertools.product([True,False],repeat=len(others_input)):
                truth_number = 0
                key=list(other_state)
                #set 0 , 0
                for inpu in first[0]:
                    index=inputs.index(inpu)
                    if inpu in positives:
                        key[index]=False
                    else:
                        key[index]= True
                for inpu in second[0]:
                    index=inputs.index(inpu)
                    if inpu in positives:
                        key[index] = False
                    else:
                        key[index] = True

                if self.logical_dict[tuple(key)]:
                    truth_number +=1
                #set 0 , 1
                for inpu in second[0]:
                    index=inputs.index(inpu)

                    if inpu in positives:
                        key[index] = True
                    else:
                        key[index] = False
                if self.logical_dict[tuple(key)]:
                    truth_number +=1
                #set 1 , 0

                for inpu in first[0]:
                    index=inputs.index(inpu)

                    if inpu in positives:
                        key[index] = True
                    else:
                        key[index] = False

                for inpu in second[0]:
                    index=inputs.index(inpu)

                    if inpu in positives:
                        key[index] = False
                    else:
                        key[index] = True
                if self.logical_dict[tuple(key)]:
                    truth_number +=1
                #set 1, 1

                for inpu in second[0]:
                    index=inputs.index(inpu)

                    if inpu in positives:
                        key[index] = True
                    else:
                        key[index] = False
                if self.logical_dict[tuple(key)]:
                    truth_number +=1

                if truth_number == 1:
                    operation_relation.append((pair[0],pair[1],'and'))
                    break
                elif truth_number == 3:
                    operation_relation.append((pair[0], pair[1], 'or'))
                    break
                elif truth_number == 2:
                    operation_relation.append((pair[0], pair[1], '^'))
                    break
        flag=True
        if len(operation_relation)==0:
            #print(classes[0])
            if classes:
                clas= classes[0]
                if len(clas[0]) == 1:
                    if clas[0][0] in positives:
                        return ' ( '+clas[0][0]+' ) ', operation_relation
                    else:
                        return ' ( not ' + clas[0][0] + ' ) ', operation_relation
                else:
                    inputs=clas[0]
                    temp = ' ( '
                    for i in range(len(inputs)):
                        input = inputs[i]

                        if input in positives:
                            if i == len(inputs) - 1:
                                temp = temp + input + ' ) '
                            else:
                                temp = temp + input + ' ' + clas[1] + ' '

                        else:
                            if i == len(inputs) - 1:
                                temp = temp + ' not ' + input + ' ) '
                            else:
                                temp = temp + ' not ' + input + ' ' + clas[1] + ' '
                return temp, operation_relation
            return '0-input',operation_relation

        flag, relation = find_final_relation(operation_relation, inputs)
        count=0
        while flag:
            operation_relation = condense_relation(operation_relation)
            flag, relation = find_final_relation(operation_relation, inputs)
            count +=1
            if count > len(classes):
                #print('Not found')
                #print(relation)
                return 'not found', operation_relation
        #print(self.positive_signals)

        function=construct_function(relation,positives)
        for state in itertools.product([True,False],repeat=len(inputs)):
            Conditions={}
            for i in range(len(inputs)):
                Conditions[inputs[i]]=state[i]
            if eval(function,Conditions) != self.logical_dict[tuple(state)]:
                return function, operation_relation
        return function,operation_relation
    def update_gene(self, Conditions):
        condition=[]
        for influencer in self.is_influenced:
            condition.append(Conditions[influencer.name])
        for externals in self.is_influenced_externals:
            condition.append(Conditions[externals])

        return self.logical_dict[tuple(condition)]

    def find_sensitivity(self):
        if not isinstance(self.logical_dict,dict) or self.logical_dict=={}:
            return 0,0
        number_of_inputs = len(list(self.logical_dict.keys())[0])
        truth= [True , False]
        real_size = number_of_inputs
        s = 0
        for i in range(0, number_of_inputs):
            count = 0
            for c in itertools.product(truth, repeat=number_of_inputs - 1):
                key = list(c)
                key.insert(i, False)
                temp = self.logical_dict[tuple(key)]
                key[i] = True
                if temp != self.logical_dict[tuple(key)]:
                    count += 1
            if count == 0:
                real_size -= 1
            s += count / (2 ** (number_of_inputs - 1))
        self.number_of_inputs=real_size
        return s, real_size
    def activity_bias(self):
        values=self.logical_dict.values()
        if len(values) ==0:
            return 0
        count=0
        for value in values:
            if value:
                count +=1
        return count/(len(values))
    def number_of_inputs(self):
        if not self.number_of_inputs:
            return len(self.is_influenced_externals)+len(self.is_influenced)
        return self.number_of_inputs

    def find_symmetries(self):
        if not isinstance(self.logical_dict,dict) or self.logical_dict == {}:
            return ['0'],[]

        number_of_inputs=len(list(self.logical_dict.keys())[0])

        eq_classes = []
        if self.is_influenced and self.is_influenced_externals:
            real_inputs = [i.name for i in list(self.is_influenced)]
            real_inputs = real_inputs + list(self.is_influenced_externals)
        else:
            real_inputs= [i for i in range(1,number_of_inputs+1)]
        if len(real_inputs)==1:
            return ['1'],[([real_inputs[0]],'I')]
        inputs= [i for i in range(0,number_of_inputs)]
        remaining=inputs.copy()
        while len(remaining) > 0:
            temp = remaining.pop(0)
            class_index = []
            eq_class = []
            operation = 'I'
            for i in range(len(remaining) - 1, -1, -1):
                sym, op= is_symmetric(remaining[i], temp, self)
                if sym:
                    eq_class.append(real_inputs[remaining[i]])
                    class_index.append(i)
                    operation = op
            for index in class_index:
                remaining.remove(remaining[index])
            eq_class.append(real_inputs[temp])

            eq_classes.append((eq_class,operation))

        distribution = [(len(i[0])) for i in eq_classes]
        distribution.sort()
        return distribution,eq_classes

    def behaviours(self):
        xor_count = 0
        and_count = 0
        or_count = 0
        count = 0
        conflicting=False
        if not isinstance(self.logical_dict,dict):
            return 0, False
        key_length = len(list(self.logical_dict.keys())[0])
        inputs= [i for i in range(0,key_length)]
        for pair in itertools.combinations(inputs, 2):
            other_genes = inputs.copy()
            other_genes.remove(pair[0])
            other_genes.remove(pair[1])
            and_behaviour = False
            or_behaviour = False
            xor_behaviour = False
            for other_state in itertools.product([True, False], repeat=len(other_genes)):
                Conditions = list(other_state)

                nt = 0
                if pair[0] < pair[1]:
                    Conditions.insert(pair[0],False)
                    Conditions.insert(pair[1],False)
                else:
                    Conditions.insert(pair[1],False)
                    Conditions.insert(pair[0],False)
                flag00 = self.logical_dict[tuple(Conditions)]

                if flag00:
                    nt += 1

                Conditions[pair[0]]=True
                Conditions[pair[1]]=False
                flag10 = self.logical_dict[tuple(Conditions)]
                if flag10:
                    nt += 1

                Conditions[pair[1]] = True
                Conditions[pair[0]] = False
                flag01 = self.logical_dict[tuple(Conditions)]
                if flag01:
                    nt += 1

                Conditions[pair[0]] = True
                Conditions[pair[1]] = True
                flag11 = self.logical_dict[tuple(Conditions)]

                if flag11:
                    nt += 1
                count += 1

                if flag10 and flag01 and not flag00 and not flag11:
                    xor_count += 1
                    xor_behaviour = True
                elif not flag10 and not flag01 and flag00 and flag11:
                    xor_count += 1
                    xor_behaviour = True
                elif nt == 1:
                    and_count += 1
                    and_behaviour = True
                elif nt== 3:
                    or_count +=1
                    or_behaviour = True

            if and_behaviour and or_behaviour:

                conflicting = True

        if count>0:
            return xor_count, conflicting
        return 0,conflicting
    def is_partially_symmetric(self):
        number_of_inputs = len(list(self.logical_dict.keys())[0])
        inputs = [i for i in range(0, number_of_inputs)]
        for pair in itertools.combinations(inputs,2):
            s, _ = is_symmetric(pair[1],pair[0],self)
            if s:
                return True
        return False
    def is_canalizing(self):
        influenced_genes = [i.name for i in self.is_influenced]
        influenced_externals = list(self.is_influenced_externals)
        inputs = influenced_genes + influenced_externals

        if len(inputs) < 2:
            return True
        Conditions = [False]*len(inputs)
        if len(inputs) > 9:
            return True
        for i in range(len(inputs)):
            inpu=inputs[i]
            others = inputs.copy()
            others.remove(inpu)
            Conditions[i] = True
            control =self.logical_dict[tuple(Conditions)]
            is_canalizing1 = True

            for other_state in itertools.product([True, False], repeat=len(others)):
                Conditions= list(other_state)
                Conditions.insert(i,True)
                temp = self.logical_dict[tuple(Conditions)]
                if (temp and not control) or (not temp and control):
                    is_canalizing1 = False

            Conditions[i] = False
            control = self.logical_dict[tuple(Conditions)]
            is_canalizing2 = True
            for other_state in itertools.product([True, False], repeat=len(others)):
                Conditions = list(other_state)
                Conditions.insert(i, False)
                temp = self.logical_dict[tuple(Conditions)]
                if (temp and not control) or (not temp and control):
                    is_canalizing2 = False
            if is_canalizing1:
                return True
            elif is_canalizing2:
                return True
        return False
    def construct_logical_dict(self):
        if isinstance(self.logical_dict,dict):
            return self.logical_dict
        Truth=[True,False]
        logical_dict={}
        number_of_inputs=len(self.is_influenced)+len(self.is_influenced_externals)
        is_influenced = list(self.is_influenced)
        is_influenced_externals = list(self.is_influenced_externals)
        for c in itertools.product(Truth,repeat=number_of_inputs):
            Conditions = {}
            for i in range(len(c)):
                if i<len(is_influenced):
                    Conditions[is_influenced[i].name]=c[i]
                else:
                    Conditions[is_influenced_externals[i-len(is_influenced)]]=c[i]

            new=eval(self.logical_dict,Conditions)
            logical_dict[tuple(c)]=new
        self.logical_dict=logical_dict
        return logical_dict
    def print_gene(self):
        print(" ")
        print(self.name),
        print('Influencing:'),
        print(self.is_influencing)
        print('Influenced by:'),
        print(self.is_influenced)
        print('Influenced externals:'),
        print(self.is_influenced_externals)
        print(self.logical_dict),
        print('Positives'),
        print(self.positive_signals)
        print('Negatives'),
        print(self.negative_signals)
    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

class Logical_Table_Network(Gene_Network):
    def __init__(self, name="",expressions=None,externals=None,Genes=[],attractors=[]):
        self.name = name
        self.expressions=expressions #txt file
        self.externals=externals # a list externals' names
        self.genes=Genes
        self.attractors=attractors # a list of lists. every element:[attractor itself(a list or a string), basin size]
        self.coherence=0
        self.sensitivity=0
        self.correspondingExternals= {} #keys: states(strings) for cyclec 1st element is the key
    def sensitivity_test_attractors(self):
        Genes=self.genes
        Externals= self.externals
        attractors = self.attractors
        correspondingExternals=self.correspondingExternals
        sensitivity_attractors = []
        ext_sensitivity_attractors= []
        length=0
        for attractor in attractors:
            if isinstance(attractor[0], str):  # that means steady state
                length +=1
                state = attractor[0]
                sensitivity_nodes = []
                ext_sensitivity_nodes = []
                ex_state = correspondingExternals[state]
                Genes, Conditions = self.set_conditions(state, ex_state)  # set attractors conditions
                for i in range(0, len(Genes)):
                    s,e = sensitivity_of_node(Genes[i],Genes,Externals,Conditions)
                    sensitivity_nodes.append(s)
                    if e:
                        ext_sensitivity_nodes.append(e)

            else:
                cycle = attractor[0]
                sensitivity_nodes = []
                ext_sensitivity_nodes = []
                ex_state = correspondingExternals[cycle[0]]
                length += len(cycle)
                for i in range(0, len(cycle)):
                    state = cycle[i]
                    if i == len(cycle) - 1:
                        control = cycle[0]
                    else:
                        control = cycle[i + 1]
                    Genes, Conditions = self.set_conditions(state, ex_state)  # set attractors conditions
                    for i in range(0, len(Genes)):
                        s,e = sensitivity_of_node(Genes[i],Genes,Externals,Conditions,nextState=control[i])
                        sensitivity_nodes.append(s)
                        if e:
                            ext_sensitivity_nodes.append(e)
            sensitivity_attractors.append(sum(sensitivity_nodes)/len(sensitivity_nodes))
            if ext_sensitivity_nodes:
                ext_sensitivity_attractors.append(sum(ext_sensitivity_nodes)/len(ext_sensitivity_nodes))

        sensitivity= sum(sensitivity_attractors)/len(sensitivity_attractors)
        self.sensitivity=sensitivity
        if len(ext_sensitivity_attractors) ==0 :
            return sensitivity, 0
        return sensitivity, sum(ext_sensitivity_attractors)/len(ext_sensitivity_attractors)
    def sensitivity_test(self):
        Genes = self.genes
        summ=0
        for gene in Genes:
            summ += gene.find_sensitivity()[0]
        return summ/len(Genes)

    def readData(self):
        f = open(self.expressions, "r")
        Genes = []
        Externals = []
        lines = f.readlines()
        f.close()
        for i in range(0, len(lines)):  # reading genes and its logical functions
            func = lines[i]
            words = func.split()
            func = optimize_texts(words)
            func = " ".join(words[2:])
            temp = Gene(words[0],function=func)
            Genes.append(temp)
        f = open(self.externals, "r")  # reading external components
        lines = f.readlines()
        f.close()
        for i in range(0, len(lines)):  # reading genes and its logical functions
            func = lines[i]
            words = func.split()
            words = optimize_texts(words)
            Externals.append(words[0])
        self.genes = Genes
        self.externals = Externals

        return Genes, Externals


def sensitivity_of_node(gene,Genes,Externals,Conditions,nextState=None):
    if not nextState:
        nextState = gene.activation
    sum = 0
    ext = 0
    nodesInfluencing = [i.name for i in list(gene.is_influenced)] + list(gene.is_influenced_externals)
    for influencer in nodesInfluencing:
        if influencer in Externals:
            Conditions = negate(influencer, Conditions)
            if gene.update_gene(Conditions) != nextState:
                ext += 1
                sum += 1
            Conditions = negate(influencer, Conditions)
        else:
            Conditions = negate(influencer, Conditions)
            if gene.update_gene(Conditions) != nextState:
                sum += 1
            Conditions = negate(influencer, Conditions)
    if ext == 0:
        return sum,None
    return sum, ext


def random_network(Network,method):
    def NumberOfInputs(Genes):
        l = []
        for gene in Genes:
            l.append((len(gene.is_influenced), len(gene.is_influenced_externals)))
        return l
    total=NumberOfInputs(Network.genes)
    Genes = Network.genes.copy()
    Externals = Network.externals.copy()
    truth = [True, False]
    if method!="symmetry dictated":
        for i in range(0, len(Genes)):
            if sum(total[i]) < 10:
                g=random_function(Genes,Externals,Genes[i],total[i],method)
            else:
                g=Genes[i]
            if g is None:
                print("No known method")
                break
            Genes[i]=g
        return Logical_Table_Network(Network.name + "_random", expressions=None, externals=Externals, Genes=Genes)
    elif method=="symmetry dictated":

        for i in range(0, len(Genes)):
            gene = Genes[i]
            gene.reset()
            rest = Genes.copy()
            inputs = set()
            for count in range(total[i][0]):
                random.shuffle(rest)
                inputs.add(rest.pop())
            ext_inputs = []
            gene.is_influenced = inputs

            for input in inputs:
                input.is_influencing.add(gene)

            if Externals:
                rest = Externals.copy()
                for count in range(total[i][1]):
                    random.shuffle(rest)
                    ext_inputs.append(rest.pop())
            gene.is_influenced_externals = set(ext_inputs)
            number_of_inputs = len(inputs) + len(ext_inputs)
            if number_of_inputs == 0:
                gene.number_of_truths = 0
                continue
            if number_of_inputs == 1:
                random.shuffle(truth)
                logical_dict = {}
                logical_dict[(True,)] = truth[0]
                logical_dict[(False,)] = truth[1]
                gene.logical_dict = logical_dict
                gene.number_of_truths = 1
                continue
            partition = choose_random_partition(number_of_inputs)
            eq_classes = []
            rest = [x for x in range(number_of_inputs)]
            for eq_class_size in partition:
                eq_class = random.sample(rest, k=eq_class_size)
                negations = random.choices([False], k=eq_class_size)
                eq_classes.append((eq_class, negations))
                rest = list(set(rest).difference(set(eq_class)))
            flag = True

            while flag:
                logical_dict = {}
                desired_number_of_truths = random.choice([i for i in range(1, 2 ** number_of_inputs)])
                all = [x for x in range(number_of_inputs)]
                count = 0
                number_of_truths = 0
                for c in itertools.product(truth, repeat=len(all)):
                    key = list(c)

                    temp = dictation(key, eq_classes, logical_dict)
                    if temp == 'new':
                        if count == 0:
                            t = random.choice(truth)
                            logical_dict[c] = t
                            count = 1
                            if t:
                                number_of_truths += 1
                        elif number_of_truths >= desired_number_of_truths:
                            logical_dict[c] = False
                            count += 1
                        else:
                            logical_dict[c] = True
                            number_of_truths += 1
                            count += 1

                    else:
                        logical_dict[c] = temp
                        count += 1
                        if temp:
                            number_of_truths += 1
                gene.logical_dict = logical_dict
                if abs(desired_number_of_truths - number_of_truths) > 1:
                    desired_number_of_truths = random.choice([i for i in range(1, 2 ** number_of_inputs)])
                s, real_size = gene.find_sensitivity()
                if real_size != number_of_inputs:
                    continue
                if method == 'xor excluded':
                    xor_count, conf = gene.behaviours()
                    if xor_count != 0 or conf != False:
                        continue
                flag = False
                gene.number_of_truths = number_of_truths
            temp = [False] * (len(inputs) + len(gene.is_influenced_externals))
            for c in range(0, len(inputs)):
                temp[c] = True

                if logical_dict[tuple(temp)] == False:
                    gene.negative_signals.add(list(inputs)[c])
                else:
                    gene.positive_signals.add(list(inputs)[c])
                temp[c] = False
        return Logical_Table_Network(Network.name + "_random", expressions=None, externals=Externals, Genes=Genes)

"""
Utilities for symmetry dictation
"""





def random_function(Genes,Externals,gene,input_size,method):
    truth = [True, False]
    gene.reset()
    rest = Genes.copy()
    inputs = set()
    for count in range(input_size[0]): #random inputs
        random.shuffle(rest)
        inputs.add(rest.pop())
    ext_inputs = []
    gene.is_influenced = inputs
    for input in inputs:
        input.is_influencing.add(gene)

    if Externals:
        rest = Externals.copy()
        for count in range(input_size[1]):
            random.shuffle(rest)
            ext_inputs.append(rest.pop())
    gene.is_influenced_externals = set(ext_inputs)
    if len(inputs) + len(ext_inputs) == 0: #if there is no input
        gene.number_of_truths = 0
        gene.logical_dict = {}
        return gene
    if len(inputs) + len(ext_inputs) == 1: # if there is one input
        random.shuffle(truth)
        logical_dict = {}
        logical_dict[(True,)] = truth[0]
        logical_dict[(False,)] = truth[1]
        gene.logical_dict = logical_dict
        gene.number_of_truths = 1
        return gene
    if method=="ncf":
        canalizing_order = list(inputs) + ext_inputs
        random.shuffle(canalizing_order)
        lines = []
        for i in range(len(canalizing_order) - 1):
            canalizing_input = canalizing_order[i]
            lines.append((canalizing_input, random.choice(truth),
                          random.choice(truth)))  # input,canalized variable, canalized value
        lines.append((canalizing_order[-1], random.choice(truth), lines[-1][2]))
        lines.append(("default", not lines[-1][2]))
        logical_dict = {}
        for c in itertools.product(truth, repeat=len(canalizing_order)):
            for i in range(len(lines)):
                line = lines[i]

                canalizing_input = line[0]
                if canalizing_input == "default":
                    logical_dict[c] = line[1]
                    break
                canalizing_value = line[1]
                canalized_value = line[2]
                if canalizing_value == c[i]:
                    logical_dict[c] = canalized_value
                    break
        gene.logical_dict = logical_dict
        return gene
    flag = True  #generate until condition is satisfied
    while flag:
        count = 0
        logical_dict = {}
        number_of_truths = random.randint(1,2**sum(input_size)-1)
        gene.number_of_truths = number_of_truths
        truths = [True] * number_of_truths + [False] * (
                2 ** (len(inputs) + len(ext_inputs)) - number_of_truths)
        random.shuffle(truths)
        for c in itertools.product(truth, repeat=len(inputs) + len(gene.is_influenced_externals)):
            logical_dict[c] = truths[count]
            count += 1
        gene.logical_dict = logical_dict
        s, real_size = gene.find_sensitivity()
        legit_methods=["xor excluded","symmetric","asymmetric","canalizing"]
        if not real_size == len(inputs) + len(ext_inputs): #regenerate
            continue
        if method=="xor excluded" and (not gene.is_symmetric() or gene.behaviours()[0] > 0):
            continue
        if method=="symmetric" and not gene.is_symmetric():
            continue
        if method=="asymmetric" and gene.is_symmetric():
            continue
        if method=="canalizing" and len(gene.is_canalizing()[0])==0:
            continue
        if method=="nested canalizing" and not gene.is_canalizing()[1]:
            continue
        elif method not in legit_methods:
            return None
        flag = False

    temp = [False] * (len(inputs) + len(gene.is_influenced_externals))
    for c in range(0, len(inputs)):
        temp[c] = True

        if logical_dict[tuple(temp)] == False:
            gene.negative_signals.add(list(inputs)[c])
        else:
            gene.positive_signals.add(list(inputs)[c])
        temp[c] = False

    return gene

def dictation(key,eq_classes,logical_dict):
    if not logical_dict:
        'new'
    for data in eq_classes:
        eq_class=data[0]
        negations=data[1]
        specific_key=[]
        i=0
        for gene in eq_class:
            temp=key[gene]
            if negations[i]:
                temp = not temp
            i += 1
            specific_key.append(temp)

        for c in itertools.permutations(specific_key,len(eq_class)):
            checking_key=key
            i=0
            for gene in eq_class:
                checking_key[gene]=c[i]
                i+=1

            if tuple(checking_key) in logical_dict:
                return logical_dict[tuple(checking_key)]

    return 'new'

def choose_random_partition(size):
    if size==1:
        return [1]
    if size==2:
        return [2]
    if size==3:
        return [1,2]
    shuffler=['o']*size+['|']*(size-2)
    random.shuffle(shuffler)
    partition=[]
    count = 0
    for i in range(len(shuffler)):
        if shuffler[i]=='o':
            count += 1
        elif shuffler[i]=='|' and count > 0:
            partition.append(count)
            count=0
    partition.sort()
    return partition

def random_symmetric_gene(number_of_inputs):
    gene= logical_gene()
    Genes=['a','b','c','d','e','f','g','h','i','j','k','l']
    gene.is_influencing = set()  # set of genes
    gene.is_influenced = set()  # set of genes
    gene.is_influenced_externals = set()  # set of strings
    gene.negative_signals = set()
    gene.positive_signals = set()
    truth=[True, False]
    partition = choose_random_partition(number_of_inputs)
    inputs = set()
    rest = Genes.copy()
    for count in range(number_of_inputs):
        random.shuffle(rest)
        inputs.add(rest.pop())
    gene.is_influenced_externals= inputs
    gene.positive_signals=inputs
    eq_classes = []
    rest = [x for x in range(number_of_inputs)]
    for eq_class_size in partition:
        eq_class = random.sample(rest, k=eq_class_size)
        negations = random.choices([False], k=eq_class_size)
        eq_classes.append((eq_class, negations))
        rest = list(set(rest).difference(set(eq_class)))
    number_of_truths = 0
    flag = True

    while flag:
        logical_dict = {}
        desired_number_of_truths = random.choice([i for i in range(1, 2 ** number_of_inputs)])
        all = [x for x in range(number_of_inputs)]
        count = 0
        number_of_truths = 0
        for c in itertools.product(truth, repeat=len(all)):
            key = list(c)

            temp = dictation(key, eq_classes, logical_dict)
            if temp == 'new':
                if count == 0:
                    t = random.choice(truth)
                    logical_dict[c] = t
                    count = 1
                    if t:
                        number_of_truths += 1
                elif number_of_truths >= desired_number_of_truths:
                    logical_dict[c] = False
                    count += 1
                else:
                    logical_dict[c] = True
                    number_of_truths += 1
                    count += 1

            else:
                logical_dict[c] = temp
                count += 1
                if temp:
                    number_of_truths += 1
        gene.logical_dict = logical_dict
        if abs(desired_number_of_truths - number_of_truths) > 1:
            desired_number_of_truths = random.choice([i for i in range(1, 2 ** number_of_inputs)])
        s, real_size = gene.find_sensitivity()
        if real_size != number_of_inputs:
            continue
        xor_count, conf = gene.behaviours()
        if xor_count != 0 or conf != False:
            continue
        flag = False
    gene.number_of_truths=number_of_truths
    return gene

def calculate_effect(gene,influencer,nodes):
    nodesInfluencing=nodes.copy()
    nodesInfluencing.remove(influencer)
    externals=list(gene.is_influenced_externals)
    if not nodesInfluencing:
        return 1
    size=len(nodesInfluencing)
    temp=[True,False]
    sum=0
    Conditions={}
    for c in itertools.product(temp,repeat=size):
        for i in range(0,size):
            if i<len(nodesInfluencing):
                node= nodesInfluencing[i]
            else:
                node= externals[i-len(nodesInfluencing)]
            Conditions[node]=c[i]
        Conditions[influencer]=False
        temp=gene.update_gene(Conditions)
        Conditions[influencer] = True
        if temp != gene.update_gene(Conditions):
            sum += 1
    return float(sum)/(2**size)


def is_symmetric(input1,input2,gene): #checks whether 2 inputs are symmetric
    number_of_inputs = len(list(gene.logical_dict.keys())[0])
    inputs = [i for i in range(0, number_of_inputs)]
    Truth=[True,False]
    inputs.remove(input1)
    inputs.remove(input2)
    input1_is_positive=False
    Falses = [False] * number_of_inputs
    Falses[input2] = True
    operation = 'None'
    sym=True
    cross_sym=True

    for c in itertools.product(Truth, repeat=len(inputs)):
        Conditions= [False]*number_of_inputs
        j=0
        for i in range(0,number_of_inputs):
            if i==input2 or i==input1:
                continue
            Conditions[i]=c[j]
            j +=1
        truth_number=0
        if operation=='None':
            Conditions[input1] = False
            Conditions[input2] = False
            if gene.logical_dict[tuple(Conditions)]:
                truth_number +=1

            Conditions[input1] = True
            Conditions[input2] = False
            if gene.logical_dict[tuple(Conditions)]:
                truth_number += 1

            Conditions[input1] = False
            Conditions[input2] = True
            if gene.logical_dict[tuple(Conditions)]:
                truth_number += 1

            Conditions[input1] = True
            Conditions[input2] = True
            if gene.logical_dict[tuple(Conditions)]:
                truth_number += 1
            if truth_number==1:
                operation='and'
            elif truth_number==2:
                operation='^'
            elif truth_number == 3:
                operation='or'

        if sym:
            Conditions[input1]=True
            Conditions[input2] = False

            temp = gene.logical_dict[tuple(Conditions)]

            Conditions[input1] = False
            Conditions[input2] = True

            if gene.logical_dict[tuple(Conditions)] != temp:
                sym = False
        if cross_sym:
            Conditions[input1] = True
            Conditions[input2] = True
            temp = gene.logical_dict[tuple(Conditions)]
            Conditions[input1] = False
            Conditions[input2] = False
            if gene.logical_dict[tuple(Conditions)] != temp:
                cross_sym=False
    t = sym or cross_sym
    return t, operation

"""
Utilities for symmetry detection and function construction
"""

def search_operation(first,second,relation):
    for r in relation:
        if r[0] == first and second == r[1]:
            return r
        if r[1] == first and second == r[0]:
            return r

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    if lst3:
        return lst3[0]
    return lst3

def condense_relation(relation):
    new_relation=relation.copy()
    flag=False
    for pair in itertools.combinations(relation,2):
        first=pair[0]
        second=pair[1]
        first_elements=[]
        second_elements=[]
        first_op=0
        second_op=0

        for data in first:
            if isinstance(data,str):
                first_op=data
            else:
                first_elements.append(data)

        for data in second:
            if isinstance(data, str):
                second_op = data
            else:
                second_elements.append(data)

        common= intersection(first_elements,second_elements)
        if not bool(common):
            continue
        elif first_op == second_op:
            first_elements.remove(common)
            second_elements.remove(common)

            r=search_operation(first_elements[0],second_elements[0],relation)
            if r:
                temp=(common,r,first_op)
                if temp not in new_relation:
                    new_relation.append(temp)
                    flag=True

    return new_relation

def find_final_relation(relation,inputs):
    inputs=set(inputs)
    for r in relation:
        data=list(r).copy()
        flag=True
        new_all=[]
        while flag:
            flag = False
            for dat in data:
                if isinstance(dat,str):
                    new_all.append(dat)
                else:
                    new_all.extend(dat)
                    flag = True
            data=new_all
            new_all=[]
        if set(data).issuperset(inputs):
            return False,r
    return True,relation

def construct_function(relation,positives):
    if isinstance(relation[0],list):
        inputs=relation[0]
        temp = ' ( '
        if len(inputs)==1:
            if inputs[0] in positives:
                return '  '+inputs[0]+ '  '
            return '  not ' + inputs[0] + '  '
        for i in range(len(inputs)):
            input=inputs[i]
            if input in positives:
                if i==len(inputs)-1:
                    temp = temp + input + ' ) '
                else:
                    temp = temp + input + ' '+relation[1]+' '

            else:
                if i==len(inputs)-1:
                    temp = temp + ' not '+ input + ' ) '
                else:
                    temp = temp + ' not '+ input + ' '+relation[1]+' '
        return temp
    else:
        return ' ( '+construct_function(relation[0],positives) + ' '+ relation[2]+' ' + construct_function(relation[1],positives)+'  )'

def negate(node,Conditions):
    if Conditions[node]:
        Conditions[node]=False
    else:
        Conditions[node] = True
    return Conditions

def to_all_string(relation):
    new_all1 = []
    flag = True
    while flag:
        flag = False
        for dat in relation:
            if isinstance(dat, str):
                new_all1.append(dat)
            else:
                new_all1.extend(dat)
                flag = True
        relation = new_all1
        new_all1 = []
    return relation

"""
Utilities for monotonicity analysis
"""

def restrict(log_dict,index,truth):
    new_log_dict={}
    Truth = [0, 1]
    key_length = int(np.log2(log_dict.size))
    for key in itertools.product(Truth, repeat=key_length):
        if key[index]==truth:
            key1=list(key)
            key1.pop(index)
            key1=tuple(key1)
            new_log_dict[key1]=log_dict[key]
    return new_log_dict


def compare_logic_dict(dict1,dict2):
    temp=[True,True]
    for key in dict1:
        if dict1[key] and not dict2[key]:    #dict1 < dict2
            temp[0]=False
            break
    for key in dict2:
        if dict2[key] and not dict1[key]:  #dict1 > dict2
            temp[1]=False
            break
    return temp



