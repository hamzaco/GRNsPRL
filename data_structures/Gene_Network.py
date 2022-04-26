from data_structures.Gene import Gene
import random
from simulation import simulate
import itertools


class Gene_Network:
    def __init__(self, name="",expressions=None,externals=None,Genes=[],attractors=[],corresponding_externals=[]):
        self.name = name
        self.expressions=expressions #txt file
        self.externals=externals # a list externals' names
        self.genes=Genes
        self.attractors=attractors # a list of lists. every element:[attractor itself(a list or a string), basin size]
        self.coherence=0
        self.sensitivity=0
        self.correspondingExternals= [] #keys: states(strings) for cyclec 1st element is the key
        if expressions:
            self.readData()
        self.initialize_effects()

    def readData(self): #reading data from expressions and externals txt file
        f = open(self.expressions, "r+")
        Genes = []
        Externals = []
        lines = f.readlines()
        f.close()
        for i in range(0, len(lines)):  # reading genes and its logical functions
            func = lines[i]
            words = func.split()
            func = optimize_texts(words)
            func = " ".join(words[2:])
            temp = Gene(words[0], False, func)
            Genes.append(temp)
        f = open(self.externals, "r")  # reading external components
        lines = f.readlines()
        f.close()
        for i in range(0, len(lines)):  # reading genes and its logical functions
            func = lines[i]
            words = func.split()
            words = optimize_texts(words)
            Externals.append(words[0])
        self.genes=Genes
        self.externals = Externals
        return Genes, Externals

    def initialize_effects(self):  # initialize network relations
        Genes=self.genes

        for gene in Genes:

            gene.is_influenced=set()
            gene.is_influenced_externals=set()
            gene.is_influencing=set()
            gene.positive_signals=set()
            gene.negative_signals=set()

        for index in range(0, len(Genes)):
            Conditions = {}
            g = Genes[index]
            words = g.function.split()
            not_counter = []
            for word in words:
                if word == 'not':
                    not_counter.append(0)
                elif word == '(':
                    if not_counter:
                        not_counter[-1] = not_counter[-1] + 1
                elif word == ')':
                    if not_counter:
                        not_counter[-1] = not_counter[-1] - 1
                        if not_counter[-1] == 0:
                            not_counter.pop()
                elif not (word == 'and' or word == 'or' or word == 'True'):
                    a = search_gene(word, Genes)
                    if a != -1:
                        Genes[a].is_influencing.add(g)
                        g.is_influenced.add(Genes[a])
                        if len(not_counter) % 2 == 1:
                            if word not in g.positive_signals:
                                g.negative_signals.add(word)
                        else:
                            if word not in g.negative_signals:
                                g.positive_signals.add(word)
                    else:
                        g.is_influenced_externals.add(word)
                        if len(not_counter) % 2 == 1:
                            if word not in g.positive_signals:
                                g.negative_signals.add(word)
                        else:
                            if word not in g.negative_signals:
                                g.positive_signals.add(word)
                    if not_counter:
                        if not_counter[-1] == 0:
                            not_counter.pop()

            g.in_degree=len(g.is_influenced_externals)+len(g.is_influenced)
        self.genes=Genes
        return Genes

    def set_conditions(self, state,exstate=None): #updates genes given a specific state(string of 0 and 1)

        Genes=self.genes
        Conditions = {}
        for i in range(0, len(state)):
            if state[i] == '1':
                Genes[i].activate()
            else:
                Genes[i].deactivate()
            Conditions[Genes[i].name] = Genes[i].activation
        self.genes=Genes
        if exstate:
            for i in range(0,len(self.externals)):
                if exstate[i]=='1':
                    Conditions[self.externals[i]]=True
                else:
                    Conditions[self.externals[i]]=False
        return Genes, Conditions

    def coherence_network(self):
        sum = 0

        if not self.attractors:
            self.calculate_attractors()
        for attractor in self.attractors:
            if isinstance(attractor[0], str):
                sum += (attractor[1] / 1000) * self.coherence_of_state(attractor[0])
            else:
                temp = 0
                for state in attractor[0]:
                    temp += self.coherence_of_state(state) / len(attractor[0])
                sum += temp * (attractor[1] / 1000)
        self.coherence=sum
        return sum

    def coherence_of_state(self,attractor):
        Genes, Conditions = self.set_conditions(attractor)
        sum = 0
        for gene in Genes:
            n = 0
            p = 0

            for negative in gene.negative_signals:
                if Conditions[negative.name]:
                    n += 1
            for positive in gene.positive_signals:
                if Conditions[positive.name]:
                    p += 1
            neg=list(gene.negative_signals)
            pos=list(gene.positive_signals)
            for i in range(0,len(gene.negative_signals)):
                if neg[i] in pos:
                    n -=1
                    p -=1
            if min(n, p) == 0:
                sum += 1
        return sum / len(Genes)

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
    def external_sensitivity(self):
        Genes=self.genes
        exts=[]
        for gene in Genes:
            extsum=0
            if len(gene.is_influenced_externals)==0:
                continue
            nodesInfluencing = [i.name for i in list(gene.is_influenced)]+list(gene.is_influenced_externals)
            for ext in gene.is_influenced_externals:
                extsum += calculate_effect(gene,ext,nodesInfluencing)
            extsum=extsum
            exts.append(extsum)
        if len(exts)==0:
            return 0
        return sum(exts)/len(exts)

    def sensitivity_test(self):
        Genes=self.genes
        Externals = self.externals
        summ = 0
        for gene in Genes:
            summ += gene.find_sensitivity()[0]
        self.sensitivity = float(summ) / len(Genes)
        return float(summ) / len(Genes)
    def calculate_attractors(self):
        [Cycles, Steady_states,correspondingExternals] = simulate(self.genes, self.externals)
        self.attractors=Cycles+Steady_states
        self.correspondingExternals=correspondingExternals
        return self.attractors

    def find_equilibrium(self,Conditions):
        States=[]
        Genes= self.genes
        eq=None
        while not eq:

            state = ""
            for i in range(0, len(Genes)):  # calculating next state
                if Genes[i].update_gene(Conditions):
                    state += "1"
                else:
                    state += "0"
            Conditions=update_conditions(Genes,Conditions)
            States.append(state)
            eq = in_equilibrium(States)
        return eq,Conditions


def random_conditions(Genes,Externals,Conditions):
    for i in range(0, len(Genes)):
        r = random.randint(0, 1)
        if r == 1:
            Conditions[Genes[i].name] = True
            Genes[i].activation = True
        else:
            Conditions[Genes[i].name] = False
            Genes[i].activation = False
    for i in range(0, len(Externals)):
        r = random.randint(0, 1)
        if r == 1:
            Conditions[Externals[i]] = True
        else:
            Conditions[Externals[i]] = False
    return Conditions



def update_conditions(Genes,Conditions):
    for gene in Genes:
        Conditions[gene.name]=gene.activation
    return Conditions



def in_steady_state(States):
    length = len(States)
    if length > 2 and States[length-1] == States[length-2]:
        #print("Steady state:")
        state = States[length-1]
        return state
    return False


def in_cycle(States):
    if len(States) > 2:
        last = States.pop()
            #print("Cycle detected:")
        if last in States:
            i = States.index(last)
            cycle = States[i:]
            #cycle.reverse()#full cycle
            return cycle
        else:
            States.append(last)
    return None


def in_equilibrium(States):
    if in_steady_state(States):
        return States[len(States)-1]
    return in_cycle(States)

def optimize_texts(text):
    for i in range(0, len(text)):
        temp = text[i]
        if temp == "AND" or temp == "OR" or temp == "NOT":
            text[i] = temp.lower()
        if ';' in list(temp):
            temp = temp.replace(";", "_")
            text[i] = temp
        if '/' in list(temp):
            temp = temp.replace("/", "_")
            text[i] = temp
        if '-' in list(temp):
            temp = temp.replace("-", "_")
            text[i] = temp
        if '+' in list(temp):
            temp = temp.replace("+","plus")
            text[i] = temp
        if '.' in list(temp):
            temp = temp.replace(".","_")
            text[i] = temp
    return text


def search_gene(gene, Genes):
    if gene == ')' or gene== '(' or gene == 'not' or gene== 'and' or gene == 'or':
        return -1
    for i in range(0, len(Genes)):
        if gene == Genes[i].name:
            return i
    return -1


def print_genes(Genes):
    for i in range(0, len(Genes)):
        Genes[i].print_gene()

def calculate_effect(gene,influencer,nodes):
    nodesInfluencing=nodes.copy()
    nodesInfluencing.remove(influencer)

    if not nodesInfluencing:
        return 1
    size=len(nodesInfluencing)
    temp=[True,False]
    sum=0
    Conditions={}

    for c in itertools.product(temp,repeat=size):
        for i in range(0,size):
            node= nodesInfluencing[i]
            Conditions[node]=c[i]
        Conditions[influencer]=False
        temp=gene.update_gene(Conditions)
        Conditions[influencer] = True
        if temp != gene.update_gene(Conditions):
            sum += 1
    return float(sum)/(2**size)

def shuffle_network(Network, method):
    Genes = Network.genes.copy()
    Externals = Network.externals.copy()
    if method == 'fixed inputs':
        for gene in Genes:
            influencers = list(gene.is_influenced.copy())
            shuffled=influencers.copy()
            random.shuffle(shuffled)
            for i in range(0, len(influencers)):
                current=influencers[i]
                gene.replace_function(current,shuffled[i])
            Genes[search_gene(gene.name,Genes)]=gene
    elif method=='swapping':
        numberOfSwaps=len(Genes)*averageNumberOfInputs(Genes)
        for i in range(0,numberOfSwaps):
            a, b = random.choices(Genes, k=2)
            while a.is_influenced.issubset(b.is_influenced) or a.is_influenced.issuperset(b.is_influenced):
                a, b = random.choices(Genes, k=2)

            temp=random.choice(list(a.is_influenced.difference(b.is_influenced)))
            temp1=random.choice(list(b.is_influenced.difference(a.is_influenced)))

            temp,temp1=a.replace_function(temp,temp1)
            temp1,temp=b.replace_function(temp1,temp)
            Genes[search_gene(a.name,Genes)]=a
            Genes[search_gene(b.name,Genes)]=b
            Genes[search_gene(temp.name,Genes)]=temp
            Genes[search_gene(temp1.name,Genes)]=temp1
    elif method == 'fixed colored swapping':
        numberOfSwaps = len(Genes) * averageNumberOfInputs(Genes)
        for i in range(0, numberOfSwaps):
            a, b = random.choices(Genes, k=2)
            positive = random.choice([True, False])
            if positive and (a.positive_signals.issubset(b.positive_signals) or a.positive_signals.issuperset(b.positive_signals)):
                continue
            if positive:
                temp = random.choice(list(a.positive_signals.difference(b.positive_signals)))
                temp1 = random.choice(list(b.positive_signals.difference(a.positive_signals)))
            if a.negative_signals.issubset(b.negative_signals) or a.negative_signals.issuperset(b.negative_signals):
                continue
            if not positive:
                temp = random.choice(list(a.negative_signals.difference(b.negative_signals)))
                temp1 = random.choice(list(b.negative_signals.difference(a.negative_signals)))

            if temp in b.is_influenced or temp1 in a.is_influenced or temp==b or temp1==a or temp==a or temp1==b or a==b:
                continue

            temp, temp1 = a.replace_function(temp, temp1)
            temp1, temp = b.replace_function(temp1, temp)
            Genes[search_gene(a.name, Genes)] = a
            Genes[search_gene(b.name, Genes)] = b
            Genes[search_gene(temp.name, Genes)] = temp
            Genes[search_gene(temp1.name, Genes)] = temp1
    elif method == 'fixed logical structure':
        for gene in Genes:
            func=gene.function
            words=func.split()
            for i in range(0,len(words)):
                if words[i] in Externals:
                    continue
                g=search_gene(words[i],list(gene.is_influenced))
                if not g==-1:
                    a=random.choice(Genes)
                    words[i]=a.name
            gene.function = " ".join(words)
        temp= Gene_Network(Network.name+"_random",None,Externals,Genes)
        temp.initialize_effects()
        return temp
    elif method== 'randomize structure':
        randomize_chance = 0.4
        for gene in Genes:
            func= gene.function
            words = func.split()
            for i in range(0,len(words)):
                word = words[i]
                r = random.random()
                if r < randomize_chance:
                    if word== 'not':
                        words[i] = ''
                    elif word == 'and':
                        words[i]='or'
                    elif word == 'or' :
                        words[i]= 'and'
                    elif word == ')' or word == '(':
                        words[i]=word
                    elif word not in Externals:
                        words[i]= 'not '+ words[i]
            gene.function= " ".join(words)

    return Gene_Network(Network.name+"_random",None,Externals,Genes)

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

def negate(node,Conditions):
    if Conditions[node]:
        Conditions[node]=False
    else:
        Conditions[node] = True
    return Conditions

def averageNumberOfInputs(Genes):
    sum = 0
    for gene in Genes:
        sum += len(gene.is_influenced)
    return sum

def shuffle_in_degrees(Network):
    for g in Network.genes:
        g.shuffle_in_degree()
