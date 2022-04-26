import time
import numpy as np

def get_partitions(collection):
    if len(collection) == 1:
        yield [ collection ]
        return

    first = collection[0]
    for smaller in get_partitions(collection[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[ first ] + subset]  + smaller[n+1:]
        # put `first` in its own subset
        yield [ [ first ] ] + smaller

class partition:
    def __init__(self,elements):
        self.elements=elements
    def __lt__(self, other):
        for element in self.elements:
            temp = False
            for o in other.elements:

                if set(o).issuperset(set(element)):
                    temp=True
            if not temp:
                return False
        return True

    def __eq__(self, other):
        for element in self.elements:
            temp = False
            for o in other.elements:
                if set(o) == set(element):
                    temp=True
            if not temp:
                return False
        return True

def A(p):
    temp=p
    if isinstance(p,partition):
        temp=p.elements
    product = 1
    for r in temp:
        product = product * (len(r) + 1)
    return 2 ** product

def supersets(p,partitions):
    summ=0
    for other in partitions:
        p1 = partition(p)
        p_other = partition(other)
        if p1 == p_other:
            continue
        if p1 < p_other:
            summ += A(p_other)
    return summ

def number_of_p_symmetric_functions(p,partitions):
    if len(p) ==1:
        return A(p)
    summ= A(p)
    for other in partitions:
        p1 = partition(p)
        p_other = partition(other)
        if p1 == p_other:
            continue
        if p1 < p_other:
            summ -= number_of_p_symmetric_functions(p_other.elements,partitions)
    return summ
Ns=np.zeros((9,9))

for n in range(3,7):
    t=time.time()
    something = list(range(1,n+1))
    partitions=get_partitions(something)
    number_of_all_symmetric_functions=0
    temp=[]
    distributions={}
    for p in partitions:
        temp.append(p)
    partitions=temp

    l=len(partitions)
    for p in partitions:
        i=0
        dist = [len(i) for i in p]
        dist.sort()
        dist = tuple(dist)
        if dist in distributions:
            multiplicity = distributions[dist][1] + 1
            distributions[dist] = (p, multiplicity)
        else:
            distributions[dist] = (p, 1)
    print("Partion (p) ", "\t".expandtabs(30), "# of function", "\t".expandtabs(30), "multiplicity")
    for dist in distributions:
        p, multiplicity = distributions[dist]
        N=number_of_p_symmetric_functions(p,partitions)
        print(len(p),"- symmetric \t".expandtabs(30),N,"\t".expandtabs(30),multiplicity)
        Ns[n,len(p)] = N*multiplicity
        if len(p) < n:
            number_of_all_symmetric_functions= number_of_all_symmetric_functions + N*multiplicity
        i += 1
    #print(number_of_all_symmetric_functions)
    all_functions= 2**(2**n)

    print("Number of symmetric "+str(n)+"-input functions: ",number_of_all_symmetric_functions)
    print("Fraction of symmetric "+str(n)+"-input functions: ",number_of_all_symmetric_functions/all_functions)
    print("--------")
