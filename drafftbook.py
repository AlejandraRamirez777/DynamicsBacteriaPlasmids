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

"""
print Birth(25,173,20,10823)

print True or True
print True or False
print False or False
print "-------"
print True and True
print True and False
print False and False

print np.linspace(1,50, 50/6.25)

print range(49)

plt.hist(np.array([7,6]))
plt.title("Distribution "+str(50)+" events")
plt.ylabel("Counts (total "+str(100)+")")
plt.xlabel("Amount of plasmids -fixated-")
plt.show()
"""

Diss = np.array([3,13,7,1,22,16,1,4,3,27,8,11,9,1,9,7,1,2,26,1,1,1,12,1,12,4,6,
6,4,5,10,5,3,20,7,4,15,1,18,5,21,15,17,8,1,1,29,11,7,1])

print Diss
print Diss.size

plt.hist(Diss,bins=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30],align="left")
plt.title("Distribution "+str(50)+" events")
plt.ylabel("Counts (total "+str(50)+")")
plt.xlabel("Amount of plasmids fixated")
plt.savefig("DisK"+str(50)+".png")
