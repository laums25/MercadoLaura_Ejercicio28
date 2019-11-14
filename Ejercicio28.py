import os
import numpy as np
import matplotlib.pyplot as plt

plt.figure(1)

data = np.loadtxt("Ejer.dat")

plt.plot(data[:401,0], data[:401,1])

plt.axis('square')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig("Ejer.png")
