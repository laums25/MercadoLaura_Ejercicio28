import os
import numpy as np
import matplotlib.pyplot as plt

plt.figure(1)

data = np.loadtxt("Ejercicio28.dat")

plt.plot(data[:,0], data[:,3])

plt.axis('square')
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig("Ejercicio28.png")
