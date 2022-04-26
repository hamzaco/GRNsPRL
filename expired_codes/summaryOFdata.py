import math
import os
from readData import read_data
from simulation import simulate
from Gene_Network import print_genes
path = "C:/Users/Lenovo/PycharmProjects/INDEP/cell_collective_networks_2018"
entries = os.listdir(path)
f=open("ALL_DATA_10000sim.txt","w+")
DATA=[]
big_entries= ["CC_CD4p_Tcell_Diff",
"CC_Differentiation_T_lymphocytes",
"CC_Glucose_Repression_Signaling_2009",
"CC_HCC1954_Breast_Cell_Line_Long-term_ErbB_Network",
"CC_IGVH_mutations",
"CC_Lymphopoiesis_Regulatory_Network",
"CC_Signaling_in_Macrophage_Activation",
"CC_Signal_Transduction_in_Fibroblasts",
"CC_Treatment_of_Castration_Resistant_Prostate_Cancer",
"CC_Yeast_Apoptosis" ]
f.write("Network             Genes,  Externals, Attractors, Basin average, Distribution \n")
for i in range(0,len(entries)):
    cell = entries[i]
    if cell not in big_entries:
        print("'"+cell+"'"+ ":")
        expressions = path +"/"+cell+"/"+"expressions.ALL.txt"
        externals = path +"/"+cell+"/"+"external_components.ALL.txt"
        [Genes, Externals] = read_data(expressions, externals)
        [Cycles, Steady_states] = simulate(Genes, Externals)
        attractors=  Cycles+Steady_states
        N = len(attractors)
        numerator = 0
        denominator = 0
        for attractor in attractors:
            bi = attractor[1]
            pi = bi / N
            numerator += math.log(pi)
            denominator += pi * math.log(pi)
        Nb=numerator/denominator
        f.write(cell+", ")
        f.write(str(len(Genes))+", "+str(len(Externals))+", "+str(len(Cycles)+len(Steady_states))+", "+str(Nb)+", ")

        for attractor in attractors:
            bi = attractor[1]
            f.write(str(bi)+"/")
        f.write("\n")


f.write("('"+cell+"' ," )
        f.write("[")
        if Cycles:
            for i in range(0,len(Cycles)-1):
                cycle=Cycles[i]
                f.write("([")
                for i in range(0,len(cycle[0])-1):
                    element=cycle[0][i]
                    f.write("'"+str(element)+"'")
                    f.write(", ")
                f.write("'"+str(cycle[0][len(cycle[0])-1])+"'")
                f.write("], "+"'"+str(cycle[1])+"'"+"), \n")

            cycle = Cycles[len(Cycles)-1]
            f.write("([")
            for i in range(0, len(cycle[0]) - 1):
                element = cycle[0][i]
                f.write("'"+str(element)+"'")
                f.write(", ")
            f.write("'"+str(cycle[0][len(cycle[0]) - 1])+"'")
            if Steady_states:
                f.write("], " +"'"+ str(cycle[1])+"'" + "), \n")
            else:
                f.write("], " +"'"+ str(cycle[1])+"'" + ")]), \n")
        if Steady_states:
            for i in range(0,len(Steady_states)-1):
                state=Steady_states[i]
                f.write("( ")
                f.write("'"+str(state[0])+"'")
                f.write(", ")
                f.write("'"+str(state[1])+"'")
                f.write(")")
                f.write(",\n")

            state = Steady_states[len(Steady_states)-1]
            f.write("( ")
            f.write("'"+str(state[0])+"'")
            f.write(", ")
            f.write("'"+str(state[1])+"'")
            f.write(")")
            f.write("]),\n")
