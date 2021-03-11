import numpy as np
import matplotlib.pyplot as plt

Avg = np.array([20.08,19.74,19.18,18.69,18.65,18.90])
StD = np.array([0.752,0.691,0.720,0.692,0.745,0.723])
La_Si = np.array([25,50,75,100,200,300])
Color = np.array(["darkviolet","blue","green","orange","red","k"])
y=np.array([1,2,3,4,5,6])

cou=0
for i in Avg:
    C = Color[cou]
    plt.errorbar(Avg[cou],y[cou],xerr=StD[cou],c=C,label= "Size:"+str(La_Si[cou]))
    plt.scatter(Avg[cou],y[cou],c=C)
    plt.legend(loc=1,fontsize="small")
    plt.ylim((-3,10))
    plt.yticks([],[])
    plt.xlabel("Average values")
    plt.title("Average values of each size with Standard Deviation")
    cou+=1

plt.savefig("Avg_Dis.png")
