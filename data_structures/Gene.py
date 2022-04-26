import itertools
import numpy as np
import random
class Gene:

    def __init__(self, name="", activation=0, function=""):
        self.name = name
        self.activation = activation #boolean
        self.function = function
        self.is_influencing = set() #set of genes
        self.is_influenced = set()  #set of genes
        self.is_influenced_externals = set() #set of strings
        self.negative_signals = set() #set of genes
        self.positive_signals= set()  # set of genes
        self.sensitivity=0
        self.dist=[]
        self.classes=[]
        self.in_degree=len(self.is_influenced_externals)+len(self.is_influenced)
        self.structure=''
        self.logical_array=None
        self.logical_dict=None
    def add_logical_function(self,f):
        self.function = f

    def activate(self):
        self.activation = True

    def deactivate(self):
        self.activation = False

    def construct_logical_array(self):
        Truth=[True,False]
        logical_array=np.ndarray(tuple([2]*(len(self.is_influenced_externals)+len(self.is_influenced))))
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
            temp=[int(key) for key in c]

            new=self.update_gene(Conditions)
            logical_array[tuple(temp)]=int(new)
        self.logical_array=logical_array
        return logical_array

    def update_gene(self, Conditions):
        if isinstance(self.logical_dict,dict):
            condition = []
            for influencer in self.is_influenced:
                condition.append(Conditions[influencer.name])
            for externals in self.is_influenced_externals:
                condition.append(Conditions[externals])

            return self.logical_dict[tuple(condition)]
        self.activation = eval(self.function, Conditions)
        return self.activation


    def find_sensitivity(self):
        if isinstance(self.logical_dict,dict):
            if not isinstance(self.logical_dict, dict) or self.logical_dict == {}:
                return 0, 0
            number_of_inputs = len(list(self.logical_dict.keys())[0])
            truth = [True, False]
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
            self.number_of_inputs = real_size
            return s, real_size
        self.construct_logical_array()
        number_of_inputs = len(self.is_influenced)+len(self.is_influenced_externals)
        truth= [1 , 0]
        real_size = number_of_inputs
        s = 0
        for i in range(0, number_of_inputs):
            count = 0
            for c in itertools.product(truth, repeat=number_of_inputs - 1):
                key = list(c)
                key.insert(i, 0)
                temp = self.logical_array[tuple(key)]
                key[i] = 1
                if temp != self.logical_array[tuple(key)]:
                    count += 1
            if count == 0:
                real_size -= 1
            s += count / (2 ** (number_of_inputs - 1))
        self.sensitivity=s
        self.effective_input=real_size
        return s, real_size

    def replace_function(self,current,new): #updates the function, replaces a new input for an input
        if current in self.positive_signals:
            self.positive_signals.remove(current)
            self.positive_signals.add(new)
        if current in self.negative_signals:
            self.negative_signals.remove(current)
            self.negative_signals.add(new)
        self.is_influenced.remove(current)
        self.is_influenced.add(new)
        words = self.function.split()
        for i in range(0,len(words)):
            if words[i] == current.name:
                words[i]=new.name
        self.function=" ".join(words)
        current.is_influencing.remove(self)
        new.is_influencing.add(self)
        return current,new

    def find_structure(self):
        if isinstance(self.logical_dict, dict):
            distribution, classes = self.find_symmetries()
            inputs = [i.name for i in list(self.is_influenced)]
            inputs = inputs + list(self.is_influenced_externals)
            positives = []
            for inpu in self.positive_signals:
                if isinstance(inpu, str):
                    positives.append(inpu)
                else:
                    positives.append(inpu.name)
            operation_relation = []

            for pair in itertools.combinations(classes, 2):
                Conditions = {}
                others = classes.copy()
                others.remove(pair[0])
                others.remove(pair[1])
                first = pair[0]
                second = pair[1]
                others_input = []
                for inpu in inputs:
                    others_input.append(inpu)

                for other_state in itertools.product([True, False], repeat=len(others_input)):
                    truth_number = 0
                    key = list(other_state)
                    # set 0 , 0
                    for inpu in first[0]:
                        index = inputs.index(inpu)
                        if inpu in positives:
                            key[index] = False
                        else:
                            key[index] = True
                    for inpu in second[0]:
                        index = inputs.index(inpu)
                        if inpu in positives:
                            key[index] = False
                        else:
                            key[index] = True

                    if self.logical_dict[tuple(key)]:
                        truth_number += 1
                    # set 0 , 1
                    for inpu in second[0]:
                        index = inputs.index(inpu)

                        if inpu in positives:
                            key[index] = True
                        else:
                            key[index] = False
                    if self.logical_dict[tuple(key)]:
                        truth_number += 1
                    # set 1 , 0

                    for inpu in first[0]:
                        index = inputs.index(inpu)

                        if inpu in positives:
                            key[index] = True
                        else:
                            key[index] = False

                    for inpu in second[0]:
                        index = inputs.index(inpu)

                        if inpu in positives:
                            key[index] = False
                        else:
                            key[index] = True
                    if self.logical_dict[tuple(key)]:
                        truth_number += 1
                    # set 1, 1

                    for inpu in second[0]:
                        index = inputs.index(inpu)

                        if inpu in positives:
                            key[index] = True
                        else:
                            key[index] = False
                    if self.logical_dict[tuple(key)]:
                        truth_number += 1

                    if truth_number == 1:
                        operation_relation.append((pair[0], pair[1], 'and'))
                        break
                    elif truth_number == 3:
                        operation_relation.append((pair[0], pair[1], 'or'))
                        break
                    elif truth_number == 2:
                        operation_relation.append((pair[0], pair[1], '^'))
                        break
            flag = True
            if len(operation_relation) == 0:
                # print(classes[0])
                if classes:
                    clas = classes[0]
                    if len(clas[0]) == 1:
                        if clas[0][0] in positives:
                            return ' ( ' + clas[0][0] + ' ) ', operation_relation
                        else:
                            return ' ( not ' + clas[0][0] + ' ) ', operation_relation
                    else:
                        inputs = clas[0]
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
                return '0-input', operation_relation

            flag, relation = find_final_relation(operation_relation, inputs)
            count = 0
            while flag:
                operation_relation = condense_relation(operation_relation)
                flag, relation = find_final_relation(operation_relation, inputs)
                count += 1
                if count > len(classes):
                    # print('Not found')
                    # print(relation)
                    return 'not found', operation_relation
            # print(self.positive_signals)

            function = construct_function(relation, positives)
            for state in itertools.product([True, False], repeat=len(inputs)):
                Conditions = {}
                for i in range(len(inputs)):
                    Conditions[inputs[i]] = state[i]
                if eval(function, Conditions) != self.logical_dict[tuple(state)]:
                    return function, operation_relation
            return function, operation_relation
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

            print(first)
            print(second)
            for other_state in itertools.product([True,False],repeat=len(inputs)):

                for i in range(len(inputs)):
                    Conditions[others_input[i]] = other_state[i]
                truth_number = 0
                #set 0 , 0
                for inpu in first[0]:
                    if inpu in positives:
                        Conditions[inpu]=False
                    else:
                        Conditions[inpu] = True


                for inpu in second[0]:
                    if inpu in positives:
                        Conditions[inpu]=False
                    else:
                        Conditions[inpu] = True

                if self.update_gene(Conditions):
                    truth_number +=1
                #set 0 , 1
                for inpu in second[0]:
                    if inpu in positives:
                        Conditions[inpu]=True
                    else:
                        Conditions[inpu] = False
                if self.update_gene(Conditions):
                    truth_number +=1
                #set 1 , 0

                for inpu in first[0]:
                    if inpu in positives:
                        Conditions[inpu] = True
                    else:
                        Conditions[inpu] = False

                for inpu in second[0]:
                    if inpu in positives:
                        Conditions[inpu] = False
                    else:
                        Conditions[inpu] = True
                if self.update_gene(Conditions):
                    truth_number +=1
                #set 1, 1

                for inpu in second[0]:
                    if inpu in positives:
                        Conditions[inpu] = True
                    else:
                        Conditions[inpu] = False
                if self.update_gene(Conditions):
                    truth_number +=1

                if truth_number == 1:
                    operation_relation.append((pair[0],pair[1],'and'))
                    break
                elif truth_number == 3:
                    operation_relation.append((pair[0], pair[1], 'or'))
                    break
        flag=True
        print(operation_relation)
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
                    for i in range(len(inputs)):
                        input = inputs[i]
                        temp=' ( '
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
            if eval(function,Conditions) != eval(self.function,Conditions):
                return 'not found', operation_relation
        self.structure=function
        return function,operation_relation


    def find_symmetries(self):
        """if isinstance(self.logical_dict,dict):
            if not isinstance(self.logical_dict, dict) or self.logical_dict == {}:
                return ['0'], []

            number_of_inputs = len(list(self.logical_dict.keys())[0])

            eq_classes = []
            if self.is_influenced and self.is_influenced_externals:
                real_inputs = [i.name for i in list(self.is_influenced)]
                real_inputs = real_inputs + list(self.is_influenced_externals)
            else:
                real_inputs = [i for i in range(1, number_of_inputs + 1)]
            if len(real_inputs) == 1:
                return ['1'], [([real_inputs[0]], 'I')]
            inputs = [i for i in range(0, number_of_inputs)]
            remaining = inputs.copy()
            while len(remaining) > 0:
                temp = remaining.pop(0)
                class_index = []
                eq_class = []
                operation = 'I'
                for i in range(len(remaining) - 1, -1, -1):
                    sym, op = is_symmetric(remaining[i], temp, self)
                    if sym:
                        eq_class.append(real_inputs[remaining[i]])
                        class_index.append(i)
                        operation = op
                for index in class_index:
                    remaining.remove(remaining[index])
                eq_class.append(real_inputs[temp])

                eq_classes.append((eq_class, operation))

            distribution = [(len(i[0])) for i in eq_classes]
            distribution.sort()
            return distribution, eq_classes"""
        self.classes=[]
        inputs = [i.name for i in list(self.is_influenced)]
        inputs = inputs + list(self.is_influenced_externals)
        eq_classes=[]
        remaining=inputs.copy()
        while len(remaining)>0:
            temp = remaining.pop(0)
            class_index=[]
            eq_class=[]
            for i in range(len(remaining)-1,-1,-1):
                symmetry,_=is_symmetric(remaining[i], temp, self)
                if symmetry:
                    eq_class.append(remaining[i])

                    class_index.append(i)

            for index in class_index:
                remaining.remove(remaining[index])
            eq_class.append(temp)
            eq_classes.append(eq_class)

        distribution=[len(i) for i in eq_classes]
        #self.dist=distribution
        #self.classes=eq_classes
        return distribution,eq_classes

    def behaviours(self):
        if isinstance(self.logical_dict,dict):
            xor_count = 0
            and_count = 0
            or_count = 0
            count = 0
            conflicting = False
            if not isinstance(self.logical_dict, dict):
                return 0, False
            key_length = len(list(self.logical_dict.keys())[0])
            inputs = [i for i in range(0, key_length)]
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
                        Conditions.insert(pair[0], False)
                        Conditions.insert(pair[1], False)
                    else:
                        Conditions.insert(pair[1], False)
                        Conditions.insert(pair[0], False)
                    flag00 = self.logical_dict[tuple(Conditions)]

                    if flag00:
                        nt += 1

                    Conditions[pair[0]] = True
                    Conditions[pair[1]] = False
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
                    elif nt == 3:
                        or_count += 1
                        or_behaviour = True

                if and_behaviour and or_behaviour:
                    conflicting = True

            if count > 0:
                return xor_count, conflicting
            return 0, conflicting
        xor_count = 0
        and_count = 0
        or_count=0
        count = 0
        conflicting = False
        influenced_genes = [i.name for i in self.is_influenced]
        influenced_externals = list(self.is_influenced_externals)
        inputs = influenced_externals + influenced_genes

        for pair in itertools.combinations(inputs, 2):
            other_genes = inputs.copy()
            other_genes.remove(pair[0])
            other_genes.remove(pair[1])
            and_behaviour = False
            or_behaviour = False
            xor_behaviour = False
            for other_state in itertools.product([True, False], repeat=len(other_genes)):
                Conditions = {}
                for i in range(0, len(other_genes)):
                    Conditions[other_genes[i]] = other_state[i]
                flag00 = False
                flag01 = False
                flag10 = False
                flag11 = False
                nt = 0
                Conditions[pair[0]] = False
                Conditions[pair[1]] = False

                flag00 = eval(self.function, Conditions)
                if flag00:
                    nt += 1
                Conditions[pair[0]] = True
                Conditions[pair[1]] = False

                flag10 = eval(self.function, Conditions)
                if flag10:
                    nt += 1
                Conditions[pair[0]] = False
                Conditions[pair[1]] = True

                flag01 = eval(self.function, Conditions)
                if flag01:
                    nt += 1
                Conditions[pair[0]] = True
                Conditions[pair[1]] = True

                flag11 = eval(self.function, Conditions)
                if flag11:
                    nt += 1
                if not (nt == 0 or nt == 4):
                    count += 1
                if flag10 and flag01 and not flag00 and not flag11:
                    xor_count += 1
                    xor_behaviour=True
                elif not flag10 and not flag01 and flag00 and flag11:
                    xor_count += 1
                    xor_behaviour=True
                elif nt == 1:
                    and_count += 1
                    and_behaviour=True
                elif nt==3 :
                    or_behaviour=True

            if and_behaviour and or_behaviour:
                conflicting = True
        return xor_count,conflicting

    def reset(self):
        self.is_influencing = set()  # set of genes
        self.is_influenced = set()  # set of genes
        self.is_influenced_externals = set()  # set of strings
        self.negative_signals = set()
        self.positive_signals = set()
    def activity_bias(self): #find number of truths in full truth table
        if isinstance(self.logical_dict,dict):
            values = self.logical_dict.values()
            if len(values) == 0:
                return 0
            count = 0
            for value in values:
                if value and value!='0':
                    count += 1
            return count / (len(values))
        number_of_inputs=len(self.is_influenced)+len(self.is_influenced_externals)

        """if self.logical_array!=None:
            return np.sum(self.logical_array)/(2**number_of_inputs)"""
        Truth=[True,False]
        is_influenced=list(self.is_influenced)
        is_influenced_externals = list(self.is_influenced_externals)
        count=0
        for c in itertools.product(Truth,repeat=number_of_inputs):
            Conditions={}
            for i in range(len(c)):
                if i<len(self.is_influenced):
                    Conditions[is_influenced[i].name]=c[i]
                else:
                    Conditions[is_influenced_externals[i-len(is_influenced)]]=c[i]
            if self.update_gene(Conditions):
                count += 1
        return count/(2**number_of_inputs)

    def is_canalizing(self):
        influenced_genes = [i.name for i in self.is_influenced]
        influenced_externals = list(self.is_influenced_externals)
        inputs_unchanged = influenced_genes + influenced_externals
        inputs=inputs_unchanged.copy()
        canalizing_nodes = []
        conditions = {}
        while True:
            cont = False
            if len(inputs)==0:
                return canalizing_nodes,True
            for inpu in inputs:
                others = inputs.copy()
                others.remove(inpu)
                for i in range(0, len(others)):
                    conditions[others[i]] = False

                conditions[inpu] = True

                control = self.update_gene(conditions)
                is_canalizing1 = True
                for other_state in itertools.product([True, False], repeat=len(others)):
                    for i in range(0, len(others)):
                        conditions[others[i]] = other_state[i]
                    temp = self.update_gene(conditions)
                    if temp==control:
                        continue
                    is_canalizing1 = False
                    break

                conditions[inpu] = False
                control = self.update_gene(conditions)
                is_canalizing2 = True

                for other_state in itertools.product([True, False], repeat=len(others)):
                    for i in range(0, len(others)):
                        conditions[others[i]] = other_state[i]
                    temp = self.update_gene(conditions)

                    if (temp and not control) or (not temp and control):
                        is_canalizing2 = False
                        break

                if is_canalizing1:
                    conditions[inpu] = True
                    cont=True
                    control = self.update_gene(conditions)
                    conditions[inpu] = False
                    canalizing_nodes.append((inpu,True,control))
                    inputs.remove(inpu)
                    break
                elif is_canalizing2:
                    conditions[inpu] = False
                    cont=True
                    control = self.update_gene(conditions)
                    conditions[inpu] = True
                    canalizing_nodes.append((inpu,False,control))
                    inputs.remove(inpu)

                    break

            if not cont:
                return canalizing_nodes, False
        return canalizing_nodes,True

    def is_symmetric(self):
        dist,_=self.find_symmetries()
        return not dist==[1]*self.in_degree

    def is_p_symmetric(self,p):
        dist,_ = self.find_symmetries()
        if len(dist) == p:
            return True
        return len(dist)

    def my_sensitivity(self):
        s=0
        number_of_inputs = len(self.is_influenced) + len(self.is_influenced_externals)
        function, op_rel= self.find_structure()
        dist,classes= self.find_symmetries()
        for clas in classes:
            adjacent=[]
            inputs_class=clas[0]
            s_c=1
            for relation in op_rel:
                first=relation[0]
                second=relation[1]

                if first==clas:
                    adjacent.append((second,relation[2]))
                elif second==clas:
                    adjacent.append((first,relation[2]))
            updated_adjacency=adjacent.copy()
            for i in range(len(adjacent)):
                for j in range(i+1,len(adjacent)):
                    candidate=adjacent[i]
                    other = adjacent[j]
                    new_all=[]
                    first=candidate[0]
                    second=other[0]
                    flag=True
                    first=to_all_string(first)
                    second=to_all_string(second)

                    first= set(first)
                    second= set(second)

                    if first.issuperset(second):
                        if other in updated_adjacency:
                            updated_adjacency.remove(other)
                    elif first.issubset(second):
                        if candidate in updated_adjacency:
                            updated_adjacency.remove(candidate)

            for i in updated_adjacency:
                flag=True
                inputs=to_all_string(i[0])

                func=construct_function(i[0],[])
                for e in inputs.copy():

                    if e=='or' or e=='and' or e=='not' or e=='I':
                        inputs.remove(e)
                number_of_1s=0
                for state in itertools.product([True, False],repeat=len(inputs)):
                    Condition={}
                    for index in range(len(inputs)):
                        Condition[inputs[index]]=state[index]
                    if eval(func, Condition):
                        number_of_1s +=1
                if i[1]=='and':
                    s_c = s_c * number_of_1s
                elif i[1]=='or':
                    s_c = s_c * (2**(len(inputs))-number_of_1s)

                s_c = s_c
            s_c = s_c * len(inputs_class)/(2**(number_of_inputs-1))
            s = s + s_c

        return s

    def remove_effect(self,influencer): #if an input is useless, remove it via this function
        words = self.function.split()
        if influencer in self.is_influenced_externals:
            self.is_influenced_externals.remove(influencer)
            for i in range(len(words)):
                if words[i] == influencer:
                    words[i] = 'True'
            if influencer in self.negative_signals:
                self.negative_signals.remove(influencer)
            else:
                self.positive_signals.remove(influencer)

        elif influencer.name in self.positive_signals:
            self.positive_signals.remove(influencer.name)
            self.is_influenced.remove(influencer)
            for i in range(len(words)):
                if words[i] == influencer.name:
                    words[i] = 'True'
        elif influencer.name in self.negative_signals:
            self.negative_signals.remove(influencer.name)
            self.is_influenced.remove(influencer)
            for i in range(len(words)):
                if words[i] == influencer.name:
                    words[i] = 'True'
        self.function=" ".join(words)

    def shuffle_in_degree(self,Genes,Externals,conservative=True):
        inputs = [i.name for i in self.is_influenced]
        inputs = inputs + list(self.is_influenced_externals)

        if conservative:
            temp= inputs.copy()
            random.shuffle(temp)
        else:
            random.shuffle(Genes)
            temp=[i.name for i in Genes[0:len(self.is_influenced)] ]
            self.is_influenced=set(Genes[0:len(self.is_influenced)])
            random.shuffle(Externals)
            self.is_influenced_externals=set(Externals[0:len(self.is_influenced_externals)])
            temp.extend(Externals[0:len(self.is_influenced_externals)])
            random.shuffle(temp)
        mapping={}
        for i in range(len(inputs)):
            mapping[inputs[i]]=temp[i]
        new_function=[]
        new_classes=[]
        if self.classes != []:
            for cls in self.classes:
                new_class=[]
                for element in cls[0]:
                    new_class.append(mapping[element])
                new_classes.append((new_class,cls[1]))
        self.classes=new_classes
        for e in self.function.split():
            if e in inputs:
                new_function.append(mapping[e])
            else:
                new_function.append(e)
        self.function=" ".join(new_function)
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

            new=eval(self.function,Conditions)
            logical_dict[tuple(c)]=new
        self.logical_dict=logical_dict
        return logical_dict
    def number_of_inputs(self):
        return len(self.is_influenced)+len(self.is_influenced_externals)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def print_gene(self):
        print(" ")
        print(self.name),
        print("Function: "),
        print(self.function),
        print('Negatives:'),
        print(self.negative_signals)
        print('Positives:'),
        print(self.positive_signals)
        print('Influencing:'),
        print(self.is_influencing)
        print('Influenced by:'),
        print(self.is_influenced)
        print('Influenced externals:'),
        print(self.is_influenced_externals)


"""
Utilities for symmetry detection and function re-construction
"""

def is_symmetric(input1,input2,gene): #checks whether 2 inputs are symmetric
    if gene.classes!=[]:
        for cls in gene.classes:
            if input1 in cls[0] and input2 in cls[0]:
                return True,cls[1]
        return False, None
    if isinstance(gene.logical_dict,dict):
        number_of_inputs = len(list(gene.logical_dict.keys())[0])
        inputs = [i for i in range(0, number_of_inputs)]
        Truth = [True, False]
        try:
            inputs.remove(input1)
            inputs.remove(input2)
        except:
            inputs_temp = [i.name for i in list(gene.is_influenced)]
            inputs_temp = inputs_temp + list(gene.is_influenced_externals)
            input2 = inputs_temp.index(input2)
            input1 = inputs_temp.index(input1)
            inputs.remove(input1)
            inputs.remove(input2)
        input1_is_positive = False

        Falses = [False] * number_of_inputs
        Falses[input2] = True
        operation = 'None'
        sym = True
        cross_sym = True

        for c in itertools.product(Truth, repeat=len(inputs)):
            Conditions = [False] * number_of_inputs
            j = 0
            for i in range(0, number_of_inputs):
                if i == input2 or i == input1:
                    continue
                Conditions[i] = c[j]
                j += 1
            truth_number = 0
            if operation == 'None':
                Conditions[input1] = False
                Conditions[input2] = False
                if gene.logical_dict[tuple(Conditions)]:
                    truth_number += 1

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
                if truth_number == 1:
                    operation = 'and'
                elif truth_number == 2:
                    operation = '^'
                elif truth_number == 3:
                    operation = 'or'

            if sym:
                Conditions[input1] = True
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
                    cross_sym = False
        t = sym or cross_sym
        return t, operation
    inputs=[i.name for i in list(gene.is_influenced)]
    inputs = inputs + list(gene.is_influenced_externals)
    Truth=[True,False]
    inputs.remove(input1)
    inputs.remove(input2)
    sym=True
    cross_sym=True

    operation='None'
    for c in itertools.product(Truth, repeat=len(inputs)):
        Conditions= {}
        for i in range(len(inputs)):
            if c[i]:
                Conditions[inputs[i]]=True
            else:
                Conditions[inputs[i]] = False
        truth_number=0

        if operation=='None':
            Conditions[input1] = False
            Conditions[input2] = False
            if gene.update_gene(Conditions):
                truth_number +=1

            Conditions[input1] = True
            Conditions[input2] = False
            if gene.update_gene(Conditions):
                truth_number += 1

            Conditions[input1] = False
            Conditions[input2] = True
            if gene.update_gene(Conditions):
                truth_number += 1

            Conditions[input1] = True
            Conditions[input2] = True
            if gene.update_gene(Conditions):
                truth_number += 1

            if truth_number==1:
                operation='and'
            elif truth_number == 3:
                operation='or'

        if sym:
            Conditions[input1]=True
            Conditions[input2] = False

            temp = gene.update_gene(Conditions)

            Conditions[input1] = False
            Conditions[input2] = True

            if gene.update_gene(Conditions) != temp:
                sym = False
        if cross_sym:
            Conditions[input1] = True
            Conditions[input2] = True
            temp = gene.update_gene(Conditions)
            Conditions[input1] = False
            Conditions[input2] = False
            if gene.update_gene(Conditions) != temp:
                cross_sym=False
    t= sym
    return t,operation



def condense_relation(relation):
    new_relation=relation.copy()
    for pair in itertools.combinations(relation,2):
        first=pair[0]
        second=pair[1]
        first_elements=[]
        second_elements=[]
        first_op=0
        second_op=0,

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

    return new_relation


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


def find_final_relation(relation,inputs):
    new_all = []
    inputs=set(inputs)
    #print(relation)
    for r in relation:
        data=list(r).copy()
        #print(data)
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
            #print(r)
            return False,r
        #print(data)
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


