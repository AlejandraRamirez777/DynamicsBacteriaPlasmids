import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def Gamma_Fit(name):
    #Import data
    info = np.genfromtxt(name+".txt",skip_header=1)
    #Obtain size of array
    size = info.size
    print "Size: "+str(size)
    #Obtain max value of array
    maxi = np.amax(info)
    print "Max Value: "+ str(maxi)
    #Obtain average of array
    Avg = np.average(info)
    print "Average: " + str(Avg)

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

    #Initial Histogram
    plt.hist(info, bins=B,color="b")
    plt.plot(0,0,label="Avg: "+str(Avg),c="b")
    plt.title("Distribution "+str(50)+" events | "+name+" repetitions")
    plt.ylabel("Counts (total "+str(info.size)+")")
    plt.xlabel("Amount of plasmids -fixated-")
    plt.legend(loc=1,fontsize="small")


    plt.savefig("Fits/"+name+".png")

Gamma_Fit("400_4")
