

from data_structures.Gene import Gene
import random
correspondingExternals={}
Conditions={}
def random_conditions(Genes,Externals):

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




def update_conditions(Genes):
    for gene in Genes:
        Conditions[gene.name]=gene.activation


def is_new_cycle(cycle, Cycles):

    for j in range(0, len(Cycles)):
        candWithOccurence=Cycles[j]
        cand= candWithOccurence[0]
        for i in range(0,len(cand)):
            temp = cand[1:]
            temp.append(cand[0])
            cand = temp
            if cand == cycle:
                Cycles.remove(candWithOccurence)
                Cycles.append((candWithOccurence[0],candWithOccurence[1]+1))
                return False
    return True


def in_steady_state(Steady_states, States,ex_state):
    length = len(States)
    if length > 2 and States[length-1] == States[length-2]:
        #print("Steady state:")
        state = States[length-1]
        for i in range(0, len(Steady_states)):
            temp = Steady_states[i]
            if temp[0] == state:
             #   print(state)
                Steady_states.remove(temp)
                Steady_states.append((temp[0], temp[1]+1))
                return True
        Steady_states.append((state, 1))
        correspondingExternals[state]=ex_state
        return True
    return False


def in_cycle(Cycles,States,ex_state):
    if len(States) > 2:
        last = States.pop()
            #print("Cycle detected:")
        if last in States:
            i = States.index(last)
            cycle = States[i:]
            #cycle.reverse()#full cycle
            if is_new_cycle(cycle, Cycles):
                #print(cycle)

                Cycles.append((cycle,1))
                correspondingExternals[cycle[0]] = ex_state
                return True
            else:
                return True
        else:
            States.append(last)
    return False


def in_equilibrium(Cycles, Steady_states, States,ex_state):
    if in_steady_state(Steady_states,States,ex_state):
        return True
    elif in_cycle(Cycles, States,ex_state):
        return True
    return False


def simulate(Genes, Externals,numIteration=0,given_state=None):
    States = []
    Cycles = []
    Steady_states = []
    numRandomVariables= len(Genes)+len(Externals)
    if numIteration==0:
        numIteration=min(pow(2, numRandomVariables+1),1000)
    if given_state:
        numIteration=1
    #Genes = iniializ,e_effects(Genes, Externals) #only for creating influence information for each gene
    for iteration in range(0, numIteration):
        if not given_state:

            random_conditions(Genes, Externals)
        ex_state = ""
        for e in Externals:
            if Conditions[e]:
                ex_state += "1"
            else:
                ex_state += "0"

        while not in_equilibrium(Cycles, Steady_states, States,ex_state):
            state = ""
            current=Conditions.copy()
            for i in range(0, len(Genes)):   #calculating next state
                if Genes[i].update_gene(current):
                    state += "1"
                else:
                    state += "0"
            update_conditions(Genes)
            States.append(state)
            if len(States) > 1000:
                break
        States = []
    return Cycles, Steady_states, correspondingExternals
