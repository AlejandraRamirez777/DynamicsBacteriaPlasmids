import numpy as np
import matplotlib.pyplot as plt
import math

#Import data
info = np.genfromtxt("DisF.txt",usecols=0,skip_header=1)

#Obtain size of array and max value
size = info.size
print "Size: "+str(size)
maxi = np.amax(info)
print "Max Value: "+ str(maxi)

#Generate the bins to clearly visualize histogram
B = np.array([])
for i in range(int(maxi)):
    B = np.append(B,i)

#Obtains the most frequent value
info_int = info.astype(int)
freq = np.bincount(info_int).argmax()
print "Most freq: "+str(freq)

plt.hist(info, bins=B)
plt.plot(0,0,label="Most freq: "+str(freq),c="b")
plt.title("Distribution "+str(50)+" events | "+str(400)+" repetitions")
plt.ylabel("Counts (total "+str(info.size)+")")
plt.xlabel("Amount of plasmids -fixated-")
plt.legend(loc=1,fontsize="small")
#plt.xlim((0,90))
plt.savefig("50eve_400rep.png")
