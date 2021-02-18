import numpy as np
import matplotlib.pyplot as plt

def Calculate():
    Col_Name = np.genfromtxt("Avgs.txt", usecols=0, dtype="string")
    Col_Info = np.genfromtxt("Avgs.txt", usecols=1)

    Info_100 = np.array([])
    Info_150 = np.array([])
    Info_200 = np.array([])
    Info_250 = np.array([])
    Info_300 = np.array([])
    Info_350 = np.array([])
    Info_400 = np.array([])

    CC = 0
    for i in Col_Name:
        if "100" in i:
            Info_100 = np.append(Info_100,Col_Info[CC])
        elif "150" in i:
            Info_150 = np.append(Info_150,Col_Info[CC])
        elif "200" in i:
            Info_200 = np.append(Info_200,Col_Info[CC])
        elif "250" in i:
            Info_250 = np.append(Info_250,Col_Info[CC])
        elif "300" in i:
            Info_300 = np.append(Info_300,Col_Info[CC])
        elif "350" in i:
            Info_350 = np.append(Info_350,Col_Info[CC])
        elif "400" in i:
            Info_400 = np.append(Info_400,Col_Info[CC])
            CC+=1

    Mega = np.array([Info_100,Info_150,Info_200,Info_250,Info_300,Info_350,
Info_400])


    for array in Mega:
        AvgL = np.array([])
        Dif = np.array([])
        for i in array:
            Avg_i = array[array != i].mean()
            AvgL = np.append(AvgL,Avg_i)

        for k in np.arange(AvgL.size):
            diff = np.absolute(AvgL[k]-array[k])
            Dif = np.append(Dif,diff)

        ans = Dif.mean()
        text_file = open("Diff.txt", "a+")
        n = text_file.write(str(ans)+" \n")
        text_file.close()

#Calculate()

def Graph():
    x = np.genfromtxt("Diff.txt", usecols=0)
    y = np.genfromtxt("Diff.txt", usecols=1)

    plt.plot(x,y)
    plt.scatter(x,y,c="k")
    plt.title("Average difference (total avg - num)")
    plt.xlabel("Repetitions Amount")
    plt.ylabel("Average difference (total avg - num)")
    plt.savefig("Differences.png")
    plt.clf()

Graph()
