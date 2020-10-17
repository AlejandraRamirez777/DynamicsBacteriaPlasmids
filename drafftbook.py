import numpy as np
import matplotlib.pyplot as plt
import math

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
