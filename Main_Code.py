import numpy as np
import matplotlib.pyplot as plt
import math
import time

#------------------------------------
#-------------FUNCTIONS--------------
#------------------------------------

#Generates Na given a Nb, utilizing the random s generated
# and thesis formula
#Param: h Hill coefficient
#Param: nnIni amount of plasmids of competition
#Return: amount of plasmids of the mutation
def NaMut(h,nnIni):
    """
    text_file = open("Graphs/Log.txt", "a+")
    n = text_file.write("nIni: "+str(nnIni)+" \n")
    text_file.close()
    """
    #Random S according to DFE
    #Return: S of according to DFE
    def randDFE():
        #Gamma distribution parameters = shape, scale
        shape=10.5
        scale=1/165.0

        ss = np.random.gamma(shape,scale)
        s= round(-1*(ss-0.040),5)
        return s

    #Turns s into Na of mutation that arised
    #Param: s Selection coefficient
    #Param: nn Nb (amount plasmids) of original
    #Param: h Hill coefficient
    #Return: Na (amount plasmids) of mutation - Integer
    def Na(s,nI,h):
        #Constant alpha of cost formula
        a = 0.029*h - 0.036
        #Na from cost formula
        Na = round(np.exp(s/a)*nI)

        #Ensure there are no bacteria without plasmids
        if Na < 0.5:
            """
            text_file = open("Graphs/Log.txt", "a+")
            n = text_file.write("Na: "+str(Na)+" \n")
            text_file.close()
            """

            R = NaMut(h,nI)

            """
            text_file = open("Graphs/Log.txt", "a+")
            n = text_file.write("R: "+str(R)+" \n")
            text_file.close()
            """
        else:
            R = Na
        return R

    SI = randDFE()
    Naa = Na(SI,nnIni,h)
    return Naa

#Calculates cost of mutation to later determine if it enters competition in bacteria level
#Param: Na of mutation
#Param: Nb of Initial
#Param: h Hill coefficient
#Return: cost of mutation and Na
def Ca(Na,Nb,h):
    #Na = Nb no enter competition because they are equal
    while (Na - Nb) == 0:
        Na = NaMut(h,Nb)

    #Constant alpha of cost formula
    a = 0.029*h - 0.036
    C = a*math.log(Nb/Na)

    return C,Na

#Defines if mutation Na dominates bacteria and enters collective competition
#Param: h Hill Coefficient
#Param: Ca cost of Mutation
#Param: Na amount of mutations
#Param: nnIni amount of original (Nb)
#Return: Na that dominates bacteria and enters collective competition
def fixed(h,CA,Na,nnIni):
    Q = np.random.random()
    #Probability fixation
    if CA >= 0:
        Prob = CA
    elif CA < 0:
        Prob = 1+CA

    #Normalization, just in case
    if Prob > 1.0:
        Prob = 1.0

    nF = 0
    if Q <= Prob:
        nF = Na

        """
        text_file = open("Graphs/Log_Mutation.txt", "a+")
        n = text_file.write("FIXED \n")
        text_file.close()
        """
    else:
        """
        text_file = open("Graphs/Log_Mutation.txt", "a+")
        n = text_file.write("NOT FIXED, TRY AGAIN \n")
        text_file.close()
        """

        Na = NaMut(h,nnIni)
        CC,Naa = Ca(Na,nnIni,h)
        nF = fixed(h,CC,Naa,nnIni)

    return nF


#--------------------------------------------------
#COMPETITION
#--------------------------------------------------

#Defines the reproduction function
#Param: NN plasmid copy number that defines bacteria
#Param: n number of bacteria with specific NN plasmids
#Return: evaluated reproduction function
def repro(NN,n):
    #Growth rate 1/25
    G = 1.0/25.0
    #Metabolic Cost plasmid (**)
    cP = 1000.0
    #Metabolic Cost bacteria
    cB = 4639221.0
    #Estimated benefit - to be measured
    b = cB/cP
    #Denominator
    D = 1 + (NN/b)
    #return
    R = (n*G)/D

    return R

#BIRTH
#Generates a reproduction decision A,B
#Param: NNa plasmid copy number that defines bacteria type A
#Param: A number of bacterias type A
#Param: NNb plasmid copy number that defines bacteria type B
#Param: B number of bacterias type B
#Return: amount of A, B after reproduction decision
#Return: probA, probB birth probabilities of decision
def Birth(NNa,A,NNb,B):

    #Normalization factor
    Pt = repro(NNa,A) + repro(NNb,B)

    #Reproduction probability type A
    probA = repro(NNa,A)/Pt

    #Reproduction probability type B
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

#DEATH
#Generates a death decision A,B
#Param: A number of bacterias type A
#Param: B number of bacterias type B
#Return: amount of A, B after death decision
def Death(A,B):

    #Normalized
    #Threshold A
    TDA = (0.5*A) / (0.5*A+0.5*B)

    #Threshold BD
    TDB = (0.5*B) / (0.5*A+0.5*B)

    #random number between 0-1
    QQ = np.random.random()

    #Either a death of A or B should happen
    #Death A
    if QQ <= TDA:
        A = A-1

    #Death B
    else:
        B = B-1

    return A,B

#Normalize
#Normalizes the amount of Bacteria to 1
#Param: AA amount of bacteria type A
#Param: BB amount of bacteria type B
#Return: cA, cB normalized amounts to 1
def nor(AA,BB):
    cA = AA / float(AA+BB)
    cB = BB / float(AA+BB)
    return cA, cB

#--------------------------------------------------
#COMPETITION
#--------------------------------------------------
#The whole process of 1 competition is performed
#Param: inA initial amount of bacterias type A
#Param: inB initial amount of bacterias type B
#Param: NNa plasmid copy number that defines bacteria type A
#Param: NNb plasmid copy number that defines bacteria type B
#Return: pA, pB arrays with data of each event for the plasmids type A and B
#Return: npA, npB normalized arrays with data of each event for the plasmids type A and B
#Return: winA, winB determines if A or B won bacteria competition (1 is win, 0 is lose)
def Go(inA, inB, NNa, NNb):
    #Stop marker to avoid infinite loops
    S = 0

    #Arrays with data of each event for the plasmids
    pA = np.array([inA])
    pB = np.array([inB])
    #Normalization
    cinA, cinB = nor(inA,inB)
    npA = np.array([cinA])
    npB = np.array([cinB])

    #Upper and lower cuts for the normalized amount of plasmids
    tUP = 0.999
    tDw = 0.001

    #Limit of number of events for the simulation
    # this limit is called LEN
    LEN = 30000

    #Markers to break loop
    C1 = True
    C2 = True

    #Code will run while the length of the array is lower than the LEN limit
    while len(npA) < LEN and S == 0:
        #Condition 1 (Upper cuts both types)
        if cinA >=tUP or cinB >= tUP:
            C1 = False

        #Condition 2 (Lower cuts both types)
        if cinA <= tDw or cinB <= tDw:
            C2 = False

        #Condition 3
        #Code will run as long as no plasmid takes completely over the population
        if (C1 == True) and (C2 == True) and (S == 0):

            #--------------BIRTH----------------
            A,B,probA,probB = Birth(NNa,inA,NNb,inB)
            inA = A
            inB = B
            pA = np.append(pA,inA)
            pB = np.append(pB,inB)

            #Normalization
            cinA, cinB = nor(inA,inB)
            npA = np.append(npA,cinA)
            npB = np.append(npB,cinB)

            #Condition 1 (Upper cuts both types)
            if (cinA >= tUP) or (cinB >= tUP):
                C1 = False
            #Condition 2 (Lower cuts both types)
            if (cinA <= tDw) or (cinB <= tDw):
                C2 = False
            #If markers to break loop have been modified, break process
            if (C1 != True) or (C2 != True):
                S = 1

            #--------------DEATH----------------
            A,B = Death(inA,inB)
            inA = A
            inB = B
            pA = np.append(pA,inA)
            pB = np.append(pB,inB)

            #Normalization
            cinA, cinB = nor(inA,inB)
            npA = np.append(npA,cinA)
            npB = np.append(npB,cinB)

            #Condition 1 (Upper cuts both types)
            if (cinA >= tUP) or (cinB >= tUP):
                C1 = False
            #Condition 2 (Lower cuts both types)
            if (cinA <= tDw) or (cinB <= tDw):
                C2 = False

            #If markers to break loop have been modified, break process
            if (C1 == False) or (C2 == False):
                S = 1

    winA = 0
    winB = 0
    #Winner count
    if len(npA) >= LEN:
        Q = np.random.random()
        probA = npA[-1]

        if Q <= probA:
            winA+=1
        else:
            winB+=1
    else:
        if npA[-1]>npB[-1]:
            winA+=1
        elif npA[-1]<npB[-1]:
            winB+=1

    #Results of the whole process
    #pA, pB arrays with amount of plasmids for each event (type A and B)
    #npA, npB arrays of pA and pB normalized to 1
    #Win of A or B (1 is win, 0 is lose)
    return pA, pB, npA, npB, winA, winB

#--------EXECUTE GO & GRAPH of COMPETITION------------
#Execute the Go function and creates the graph of the competition
#Param: inA initial amount of bacteria type A
#Param: inB initial amount of bacteria type B
#Param: NNa copy number that defines bacteria type A
#Param: NNb copy number that defines bacteria type B
#Param: IT counter distribution
#Return: winner of competition (value of NNa or NNb)
#Return: graph of the competition
def proGraph(inA, inB, NNa, NNb, IT):
    #Start measuring simulation time
    t0 = time.time()

    #Return of Go function (check Go comments)
    pA, pB, npA, npB, winA, winB = Go(inA, inB, NNa, NNb)

    """
    text_file = open("Graphs/Log.txt", "a+")
    n = text_file.write("pA: "+str(pA)+" \n")
    text_file.close()

    text_file = open("Graphs/Log.txt", "a+")
    n = text_file.write("pB: "+str(pB)+" \n")
    text_file.close()

    text_file = open("Graphs/Log.txt", "a+")
    n = text_file.write("npA: "+str(npA)+" \n")
    text_file.close()

    text_file = open("Graphs/Log.txt", "a+")
    n = text_file.write("npB: "+str(npB)+" \n")
    text_file.close()

    #Normalized graph of competition
    # type A = blue   type B = green
    plt.plot(np.linspace(0,len(npA), num = len(npA)), npA, c = "b")
    plt.plot(np.linspace(0,len(npB), num = len(npB)), npB, c = "g")

    #Not normalized graph of competition
    #plt.plot(np.linspace(0,len(pA), num = len(pA)), pA, c = "b")
    #plt.plot(np.linspace(0,len(pB), num = len(pB)), pB, c = "g")

    tsim = round(time.time()-t0,3)
    #Labels
    plt.plot(0,0, c = "b", label = "$N_{A}$ = " + str(NNa))
    plt.plot(0,0, c = "g", label = "$N_{B}$ = " + str(NNb))
    plt.plot(0,0, c = "b", linestyle = "--",label = "WinA = " + str(winA))
    plt.plot(0,0, c = "g", linestyle = "--",label = "WinB = " + str(winB))
    plt.plot(0,0, c = "y", label = "$t_{simu}$ = " + str(round(tsim,1))+"s")
    #General setting of graph
    #plt.xlim((0,400))
    #plt.ylim((0,4000))
    plt.title("Single competition $NN_{A}$="+str(NNa)+" vs $NN_{B}$="+str(NNb)+" ($n_{A}$="+str(inA)+" | $n_{B}$="+str(inB)+")",fontsize="small")
    plt.xlabel("Events")
    plt.ylabel("Normalized amount of Bacteria")
    #plt.ylabel("Amount of Bacteria")
    plt.legend(loc = 2, fontsize = "x-small",framealpha=0.8)
    plt.savefig("Graphs"+str(IT)+"/Graph_" + str(NNa) + "_" + str(NNb) + "_" + str(time.time())+ ".png")
    plt.clf()

    text_file = open("Graphs"+str(IT)+"/Results.txt", "a+")
    n1 = text_file.write("Total Win A ("+str(NNa)+") = "+str(winA)+"\n")
    n = text_file.write("Total Win B ("+str(NNb)+") = "+str(winB)+"\n")
    text_file.close()

    text_file = open("Graphs/Results"+str(IT)+".txt", "a+")
    n1 = text_file.write("Total Win A ("+str(NNa)+") = "+str(winA)+"\n")
    n = text_file.write("Total Win B ("+str(NNb)+") = "+str(winB)+"\n")
    text_file.close()
    """

    NF = -1
    if winA > winB:
        NF = NNa
    elif winA < winB:
        NF = NNb
    return NF

#--------------------------------------------------
#--------------------EXECUTE-----------------------
#--------------------------------------------------

#Execute full process with repetitions
#Param: h Hill coefficient
#Param: nnIni Initial amount of plasmids in bacteria (nB)
#Param: INN Num of bacterias that are NO mutation (NNb)
#Param: inM Num of bacterias with mutation (nA)
#Param: rep Times the competition is allowed +1 (dot graph)
#Param: disRep Times the whole process is repeated
#return: final 10 values of fixed plasmids of all disRep in bacteria
#return: graph with competitions - dot graphs
def exe(h,nnIni,INN,inM,rep,disRep):
    #Counter
    IT = 1
    #Array wins finale
    WFin = np.array([])
    #Loop for distribution repetitions
    for i in range(disRep):
        #1 INITIALIZATION ROUND

        #Mutation generation
        #Fix is the plasmid copy number of mutation
        Na = NaMut(h,nnIni)
        CC,Naa = Ca(Na,nnIni,h)
        fix = fixed(h,CC,Naa,nnIni)

        """
        text_file = open("Graphs"+str(IT)+"/Log_Mutation.txt", "a+")
        n1 = text_file.write("FIXED --> "+ str(fix)+"\n")
        n = text_file.write("----------------------------- \n")
        text_file.close()

        text_file = open("Graphs/Results"+str(IT)+".txt", "a+")
        n1 = text_file.write("Ini: "+str(nnIni)+"\n")
        n = text_file.write("Fixed: "+str(fix)+"\n")
        text_file.close()
        """

        #proGraph(inA, inB, NNa, NNb)
        #Therefore, A is the mutation
        WW = proGraph(inM, INN, fix, nnIni, IT)

        """
        text_file = open("Graphs/Results"+str(IT)+".txt", "a+")
        n1 = text_file.write("Win: "+str(WW)+"\n")
        n = text_file.write( "------------------- \n")
        text_file.close()
        """

        #Array of original NN (not the mutation)
        IniR = np.array([nnIni])
        #Mutation that dominates bacteria and enters competition
        FIX = np.array([fix])
        #Bacteria type (defined by its NN) that wins competition
        WIN = np.array([WW])

        #The one who wins starts as original in the next round
        #   to generate mutation
        nnIni = WW

        for i in range(rep):
            #Mutation generation
            #Fix is the plasmid copy number of mutation
            Na = NaMut(h,nnIni)
            CC,Naa = Ca(Na,nnIni,h)
            fix = fixed(h,CC,Naa,nnIni)

            """
            text_file = open("Graphs"+str(IT)+"/Log_Mutation.txt", "a+")
            n1 = text_file.write("FIXED --> "+ str(fix)+"\n")
            n = text_file.write("----------------------------- \n")
            text_file.close()


            text_file = open("Graphs/Results"+str(IT)+".txt", "a+")
            n1 = text_file.write("Ini: "+str(nnIni)+"\n")
            n = text_file.write("Fixed: "+str(fix)+"\n")
            text_file.close()
            """

            #proGraph(inA, inB, NNa, NNb)
            #Therefore, A is the mutation
            WW = proGraph(inM, INN, fix, nnIni, IT)

            """
            text_file = open("Graphs/Results"+str(IT)+".txt", "a+")
            n1 = text_file.write("Win: "+str(WW)+"\n")
            n = text_file.write( "------------------- \n")
            text_file.close()
            """
            #Array of original NN (not the mutation)
            IniR = np.append(IniR, nnIni)
            #Mutation that dominates bacteria and enters competition
            FIX = np.append(FIX, fix)
            #Bacteria type (defined by its NN) that wins competition
            WIN = np.append(WIN,WW)

            #The one who wins starts as original in the next round
            #   to generate mutation
            nnIni = WW

        """
        #Graphs
        plt.scatter(np.linspace(1,len(IniR), num = len(IniR)),IniR,c="k")
        plt.scatter(np.linspace(1,len(FIX), num = len(FIX)),FIX,c="k",label = "Contender")
        plt.scatter(np.linspace(1,len(WIN), num = len(WIN)),WIN,c="r",edgecolors="none",s=70,label = "Winner")
        plt.title("General competitions (Size "+str(INN+inM)+")"")
        plt.xlabel("Events")
        plt.ylabel("Plasmid copy number of bacteria")
        plt.xlim(0,len(WIN)+1)
        plt.ylim(0,np.amax(FIX)+1)
        #plt.xticks(np.linspace(1,len(WIN), num = len(WIN)/6.25))
        plt.legend(loc = 2, fontsize = "x-small",framealpha=0.8)
        plt.savefig("Graphs/Final_" + str(WIN[-1]) + "_" + str(INN) + "_"+str(IT)+".png")
        plt.clf()
        """
        IT+=1

        #-----------------------------------------
        #WARNING CHANGING VALUES
        #-----------------------------------------

        #CHANGE IF EVENTS OF DOT GRAPHS ARE CHANGED!!!
        #   (last 10 value of 50 dot graphs events)
        WFin = np.append(WFin,WIN[40:50])
        #CHANGE IF INITIAL AMOUNT OF ORIGINAL BACTERIA TYPE IS CHANGED!!!
        nnIni = 20
    return WFin

#---------------------------------------------------
#VALUES TO RUN THE SIMULATION
#---------------------------------------------------
EVE = 49
DISREP = 350
NOR = 90
MUT = 10

#exe(h,nnIni,INN,inM,rep,disRep)
#rep+1 are the number of events in the dot graph
#disRep is the number of dot graphs generated
WinDis = exe(3.0,20,NOR,MUT,EVE,DISREP)

#Print info of distribution in a DisF.txt
#   The first line is header - do not take into account in calculations
text_file = open("DisF.txt", "a+")
n1 = text_file.write("WIN FINALE DIS\n")
for i in WinDis:
    n = text_file.write(str(i)+"\n")
text_file.close()
