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
    if group == 50:
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
    plt.title("Distribution "+str(50)+" events | "+name+" repetitions")
    plt.ylabel("Counts (total "+str(info.size)+")")
    plt.xlabel("Amount of plasmids -fixated-")
    plt.legend(loc=1,fontsize="small")
    plt.savefig("Size_Fits/"+name+".png")
    plt.clf()


for i in [50,100,200]:
    for u in range(1,21):
        name = str(u)+"_S"+str(i)
        Gamma_Fit(name,i)


def Info_ana(name):
    #Import data
    info = np.genfromtxt("Size_Info/"+name+".txt",skip_header=1)

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

for i in [50,100,200]:
    for u in ["Avg","StD","Var","ErrM"]:
        name = u+"_"+str(i)
        Info_ana(name)
