import numpy as np
import matplotlib.pyplot as plt

# Number of groups
NN = 1000
#Hill coefficient
h = 3.0
#Initial n value
nnIni = 20

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
#print SI

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


#print Na(SI,nn,h)

#Defines if mutation Na gets fixed
def fixed(s,nI,h):
    Q = np.random.random()
    #Probability fixation Boettiger
    Prob = abs(2*s)
    #Normalization
    if Prob > 1.0:
        Prob = 1.0
    if Q <= Prob:
        nF = Na(s,nI,h)
    else:
        nF = nI
    return nF

nII = nnIni
for i in range(10000):
    nn = fixed(SI,nII,h)
    #print nn
    nII = nn
    #print "NN"
    #print nII
    #print "-------"

print nII
