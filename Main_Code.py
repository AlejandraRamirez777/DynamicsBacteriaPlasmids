import numpy as np
import matplotlib.pyplot as plt
import math

# Number of groups
NN = 1000
#Hill coefficient
h = 3.0
#Initial n value
nnIni = 20

#Generates probability value to determine if mutation dominates
#   bacteria and will enter competition
def CCmut(h,nnIni):
    #Random S according to DFE
    #Return: S of according to DFE
    def randDFE():
        #Gamma distribution parameters = shape, scale
        shape=10.5
        scale=1/165.0

        ss = np.random.gamma(shape,scale)
        s= round(-1*(ss-0.040),5)
        return s

    SI = randDFE()
    print "DFE_value: " + str(SI)

    #Na of mutation that arised
    #Param: s Selection coefficient
    #Param: nn Nb (amount plasmids) of competition
    #Param: h Hill coefficient
    #Return: Na (amount plasmids) of mutation - Integer
    def Na(s,nI,h):
        #Constant alpha of cost formula
        a = 0.029*h - 0.036
        #Na from cost formula
        Na = round(np.exp(s/a)*nI)
        #Ensure there are no bacteria without plasmids
        if Na < 1.0:
            Na = 1.0
        return Na

    Naa = Na(SI,nnIni,h)
    print "Na_Equivalent: "+ str(Naa)

    #Cost of mutation to determine if dominates bacteria
    #   to determine if it enters competition
    #Param: Na of mutation
    #Param: Nb of Initial
    #Param: h Hill coefficient
    #Return: cost of mutation
    def Ca(Na,Nb,h):
        #Constant alpha of cost formula
        a = 0.029*h - 0.036
        C = a*math.log(Nb/Na)
        return C

    CC = Ca(Naa,nnIni,h)
    print "CostA: "+ str(CC)

    return CC,Naa

#Defines if mutation Na dominates bacteria and enters competition
def fixed(CA,Na):
    Q = np.random.random()
    #Probability fixation
    if CA >= 0:
        Prob = CA
    elif CA < 0:
        Prob = 1+CA
    #Normalization
    if Prob > 1.0:
        Prob = 1.0

    nF = 0
    if Q <= Prob:
        nF = Na
        print "FIXED"
    else:
        print "NOT FIXED, TRY AGAIN"
        CC, Naa = CCmut(h,nnIni)
        nF = fixed(CC,Naa)
    return nF

CC, Naa = CCmut(h,nnIni)
print fixed(CC,Naa)
