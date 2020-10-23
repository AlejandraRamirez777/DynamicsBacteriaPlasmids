import numpy as np
import matplotlib.pyplot as plt
import math

"""
#Inin
A = np.array([20,23,23])
#FIXED
F = np.array([23,27,25])
#Win
WW = np.array([23,23,25])

plt.scatter(np.linspace(1,len(A), num = len(A)),A,c="k")
plt.scatter(np.linspace(1,len(F), num = len(F)),F,c="k")
plt.scatter(np.linspace(1,len(WW), num = len(WW)),WW,c="r",edgecolors="none",s=70)
plt.xlim(0,len(WW)+1)
plt.xticks(np.linspace(1,len(A), num = len(A)))
plt.show()
"""

#Defines the reproduction function
#Param: NN copy number that defines bacteria
#Param: n number of bacteria
#Return: evaluate reproduction function
def repro(NN,n):
    #Growth rate 1/25
    G = 1.0/25.0
    #Metabolic Cost plasmid (**)
    cP = 10000.0
    #Metabolic Cost bacteria
    cB = 4639221.0
    #Estimated benefit - to be measured
    b = cP/cB
    #Denominator
    D = 1 + (NN/b)
    #return
    R = (n*G)/D

    return R

#BIRTH
#Generates a reproduction decision A,B or None
#Param: NNa copy number that defines bacteria type A
#Param: A number of bacterias type A
#Param: NNb copy number that defines bacteria type B
#Param: B number of bacterias type B
#Return: amount of A, B after reproduction decision
def Birth(NNa,A,NNb,B):

    #Normalization factor
    Pt = repro(NNa,A) + repro(NNb,B)

    #Reproduction probablity type A
    probA = repro(NNa,A)/Pt

    #Reproduction probablity type B
    probB = repro(NNb,B)/Pt

    #Random number between 0-1
    Q1 = np.random.random()

    #Either a Birth of A or B should happen
    #Birth A
    if Q1 <= probA:
        A+= 1
    #Birth B
    else:
        B += 1

    return A,B,probA,probB

print Birth(25,173,20,10823)

print True or True
print True or False
print False or False
print "-------"
print True and True
print True and False
print False and False
