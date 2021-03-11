import numpy as np
import math
import scipy.stats as stats
import matplotlib.pyplot as plt

def Gamma_Fit(name,group):
    #Import data
    info = np.genfromtxt("Size_Data/"+name+".txt",skip_header=1)
    #Size of array
    size = info.size
    print "Size: "+str(size)
    #Max value of array
    maxi = np.amax(info)
    print "Max Value: "+ str(maxi)
    #Average of array
    Avg = np.average(info)
    print "Average: " + str(Avg)
    #Variance of array
    Var = np.var(info)
    #Standard deviation
    StD = math.sqrt(Var)
    #Error of the mean
    ErrM = StD/float(math.sqrt(size))

    #Print Average
    if group == 25:
        text_file = open("Size_Info/Avg_25.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_25.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_25.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_25.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S25/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()


    elif group == 50:
        text_file = open("Size_Info/Avg_50.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_50.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_50.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_50.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S50/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()

    elif group == 75:
        text_file = open("Size_Info/Avg_75.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_75.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_75.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_75.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S75/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()

    elif group == 100:
        text_file = open("Size_Info/Avg_100.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_100.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_100.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_100.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S100/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()

    elif group == 200:
        text_file = open("Size_Info/Avg_200.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_200.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_200.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_200.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S200/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()

    elif group == 300:
        text_file = open("Size_Info/Avg_300.txt", "a+")
        n = text_file.write(str(Avg)+" \n")
        text_file.close()

        text_file = open("Size_Info/Var_300.txt", "a+")
        n = text_file.write(str(Var)+" \n")
        text_file.close()

        text_file = open("Size_Info/StD_300.txt", "a+")
        n = text_file.write(str(StD)+" \n")
        text_file.close()

        text_file = open("Size_Info/ErrM_300.txt", "a+")
        n = text_file.write(str(ErrM)+" \n")
        text_file.close()

        for i in range(1,int(maxi)+1):
            count = np.count_nonzero(info==i)
            text_file = open("Size_Info/S300/Info_"+str(i)+".txt", "a+")
            n = text_file.write(str(count)+" \n")
            text_file.close()


    #----------FIT GAMMA------------------------------------------
    fit_alpha, fit_loc, fit_beta=stats.gamma.fit(info)

    print fit_alpha, fit_loc, fit_beta

    g = stats.gamma(fit_alpha, fit_loc, fit_beta)
    x = np.linspace(0,maxi,100)

    plt.plot(x,g.pdf(x)*size,c="r",lw=3,label="Alpha: "+str(round(fit_alpha,1))+
    " | Loc: "+str(round(fit_loc,1))+" | Beta: "+str(round(fit_beta,1)))

    #-------------HISTOGRAM PLOTTING------------------------
    #Generate the bins to clearly visualize histogram
    B = np.array([])
    for i in range(int(maxi)):
        B = np.append(B,i)

    #Histogram
    plt.hist(info, bins=B,color="b")
    plt.plot(0,0,label="Avg: "+str(Avg),c="b")
    plt.title("Distribution "+str(50)+" events | 350 repetitions | Size "+str(group),fontsize="small")
    plt.ylabel("Counts (total "+str(info.size)+")")
    plt.xlabel("Amount of plasmids -fixated-")
    plt.legend(loc=1,fontsize="small")
    plt.savefig("Size_Fits/"+name+".png")
    plt.clf()


#Run instruction Gamma_Fit
for i in [25,50,75,100,200,300]:
    for u in range(1,41):
        name = str(u)+"_S"+str(i)
        Gamma_Fit(name,i)


def Info_ana(name):
    #Import data
    info = np.genfromtxt("Size_Info/"+name+".txt")

    #Average of array
    Avg = np.average(info)
    #Variance of array
    Var = np.var(info)
    #Standard deviation
    StD = math.sqrt(Var)

    text_file = open("Size_Info/"+name+".txt", "a+")
    n = text_file.write("Avg_Full: "+str(Avg)+" \n")
    n2 = text_file.write("Variance: "+str(Var)+" \n")
    n3 = text_file.write("Stand_Dev: "+str(StD)+" \n")
    text_file.close()


#Run instruction Info_ana
for i in [25,50,75,100,200,300]:
    for u in ["Avg","StD","Var","ErrM"]:
        name = u+"_"+str(i)
        Info_ana(name)


def Info_Single(cc,group):
    #Import data
    info = np.genfromtxt("Size_Info/S"+str(group)+"/Info_"+str(cc)+".txt")

    #Size of array
    size1 = info.size

    print cc, group, size1

    #Fill with zeros
    Miss = 40-size1
    if Miss != 0:
        text_file = open("Size_Info/S"+str(group)+"/Info_"+str(cc)+".txt", "a+")
        for i in range(Miss):
            n = text_file.write(str(0)+" \n")
        text_file.close()

    #Reread data#Import data
    info = np.genfromtxt("Size_Info/S"+str(group)+"/Info_"+str(cc)+".txt")
    #Size of array
    size2 = info.size

    #Average of array
    Avg = np.sum(info)/float(size2)

    #Variance of array
    Var = np.var(info)
    #Standard deviation
    StD = math.sqrt(Var)

    text_file = open("Size_Info/S"+str(group)+"/Info_"+str(cc)+".txt", "a+")
    n = text_file.write(str(Avg)+" \n")
    n3 = text_file.write(str(StD)+" \n")
    text_file.close()


#Run instruction Info_Single
cA = 0
for i in [25,50,75,100,200,300]:
    endArr = np.array([191,137,187,121,134,126])
    rr = endArr[cA]
    for u in range(1,rr+1):
        Info_Single(u,i)
    cA+=1

def Grab_Values(cc,group):
    #Import data
    info = np.genfromtxt("Size_Info/S"+str(group)+"/Info_"+str(cc)+".txt")

    #Last last one is Avg
    #Last is StD
    Avg = info[-2]
    StD = info[-1]
    #Return values
    return Avg, StD

def Graph_Three(UpErr,RealV,DoErr,xR,group,StdErr,Full_Avg,cou):

    #X axis values
    x = np.linspace(1,xR,xR)

    plt.errorbar(x,RealV,yerr=StdErr,c="bisque")
    plt.plot(x,RealV,c="k",label = "Average")
    plt.plot(x,UpErr, c="orange",label = "Standard Deviation")
    plt.plot(x,DoErr, c="orange")
    plt.plot(0,0,c="white",label="Average: "+str(Full_Avg))

    #General plot
    Color = np.array(["darkviolet","blue","green","orange","red","k"])
    C = Color[cou]
    #plt.plot(x,RealV,c=C,label = "Size: "+str(group))


    plt.legend(loc=1,fontsize="small")
    plt.xlim((0,100))
    plt.ylim((0,200))
    plt.xlabel("Amount of plasmids -fixated-")
    plt.ylabel("Average counts")
    plt.title("Average distribution | 50 events | 350 repetitions | Size "+str(group), fontsize="small")
    plt.xticks([0,10,20,30,40,50,100,150,200],[0,10,20,30,40,50,100,150,200])
    #plt.xticks([0,10,20,30,40,50,100],[0,10,20,30,40,50,100])

    plt.savefig("Final_Results/"+str(group)+"_Finale.png")
    #plt.savefig("Final_Results/All_Finale_Zoom.png")
    plt.clf()

#Run instruction Graph_Final_Results
cou = 0
for i in [25,50,75,100,200,300]:
    #Create arrays for graphs
    UpErr = np.array([])
    RealV = np.array([])
    DoErr = np.array([])
    StdErr = np.array([])

    xRarr = np.array([191,137,187,121,134,126])
    xR = xRarr[cou]

    for u in range(1,xR+1):
        #Obtain values from files
        AvgI, StDI = Grab_Values(u,i)
        #Calculate error up and down
        UpV = AvgI+StDI
        DoV = AvgI-StDI
        if DoV < 0.0:
            DoV = 0.0
        #Append values to arrays
        UpErr = np.append(UpErr,UpV)
        RealV = np.append(RealV,AvgI)
        DoErr = np.append(DoErr,DoV)
        StdErr = np.append(StdErr,StDI)

    Full_Avg_Arr = np.array([20.081,19.744,19.180,18.688,18.648,18.899])
    Full_Avg = Full_Avg_Arr[cou]

    #Graph results
    Graph_Three(UpErr,RealV,DoErr,xR,i,StdErr,Full_Avg,cou)
    cou+=1
