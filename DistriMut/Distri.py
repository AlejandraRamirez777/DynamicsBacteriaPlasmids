import matplotlib.pyplot as plt
import numpy as np

#Gamma distribution parameters = shape, scale
#size --> number of outcome (for one number, leave blank)
shape=10.5
scale=1/165.0
size = 10000

ss = np.random.gamma(shape,scale,size)
s=-1*(ss-0.040)
plt.hist(s,50, normed = 0)

plt.title("DFE")
plt.ylabel("Counts (over 10,000)")
plt.xlabel("Selection Coefficient")

#Medidas
Mean = round(np.mean(s),3)
Max = round(np.max(s),3)
Min = round(np.min(s),3)
print "Mean = " + str(Mean)
print "Max = " + str(Max)
print "Min = " + str(Min)

#Legends
plt.plot(0,0,c="k",label = "Mean = "+str(Mean))
plt.plot(0,0,c="k",label = "Max = "+str(Max))
plt.plot(0,0,c="k",label = "Min = "+str(Min))
plt.vlines(0.01,0.01,800, color="magenta", linewidth=2)
plt.vlines(-0.01,0.01,800, color="magenta", linewidth=2)
plt.hlines(100, -0.149,0.049, color="lime", linewidth=2)
plt.legend(loc = 2, fontsize="small")

plt.savefig("Dis.png")
