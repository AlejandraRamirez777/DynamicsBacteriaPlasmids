import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack


def Smooth(name):
    #Import data
    info = np.genfromtxt(name+".txt",skip_header=1)

    #Obtain size of array and max value
    size = info.size
    print "Size: "+str(size)
    maxi = np.amax(info)
    print "Max Value: "+ str(maxi)
    Avg = np.average(info)
    print "Average: " + str(Avg)

    #text_file = open("Avgs.txt", "a+")
    #n = text_file.write(name +" " +str(Avg)+" \n")
    #text_file.close()

    #Generate the bins to clearly visualize histogram
    B = np.array([])
    for i in range(int(maxi)):
        B = np.append(B,i)

    #Initial Histogram
    plt.hist(info, bins=B)
    plt.plot(0,0,label="Avg: "+str(Avg),c="b")
    plt.title("Distribution "+str(50)+" events | "+name+" repetitions")
    plt.ylabel("Counts (total "+str(info.size)+")")
    plt.xlabel("Amount of plasmids -fixated-")
    plt.legend(loc=1,fontsize="small")
    #plt.xlim((0,90))
    plt.savefig("Img/Inicial_"+name+".png")
    plt.clf()

    #Fast Fourier Transform
    FT = fftpack.fft(info)
    #print FT
    sFT = FT.size
    print "Size FT: "+str(sFT)

    #Amplitude
    Amplitude = np.abs(FT)
    #Frequencies FFT
    #d defines the order of magnitude of frequency x axis
    freq = fftpack.fftfreq(sFT, d=0.001)
    print freq.size
    #print freq

    #plt.plot(fftpack.fftshift(freq),FT)
    plt.plot(freq,FT)
    plt.ylim((0,2000))
    plt.title("Fourier Transform")
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude")
    plt.savefig("Img/FT_Ini_"+name+".png")
    plt.clf()

    Amp_freq = np.array([Amplitude, freq])

    Filter = np.array([])
    for i in np.arange(sFT):
        if np.abs(freq[i]) > 300:
            Filter = np.append(Filter,0)
        else:
            Filter = np.append(Filter,Amplitude[i])

    text_file = open("Freq.txt", "a+")
    for i in freq:
        n = text_file.write(str(i)+"\n")
    text_file.close()
    print Filter
    plt.plot(freq,Filter)
    plt.ylim((0,2000))
    plt.title("Fourier Transform")
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude")
    plt.savefig("Img/Filter_FT_"+name+".png")
    plt.clf()

    #print Filter
    OK = fftpack.ifft(Filter)
    #print OK

    AvgOK = np.average(OK).real
    print "AverageOK: " + str(AvgOK)

    #Initial Histogram
    plt.hist(OK, bins=B)
    plt.plot(0,0,label="Avg: "+str(AvgOK),c="b")
    plt.title("Distribution "+str(50)+" events | "+name+" repetitions")
    plt.ylabel("Counts (total "+str(info.size)+")")
    plt.xlabel("Amount of plasmids -fixated-")
    plt.legend(loc=1,fontsize="small")
    plt.xlim((0,80))
    plt.savefig("Img/FilterOk_"+name+".png")
    plt.clf()

Smooth("400_4")
