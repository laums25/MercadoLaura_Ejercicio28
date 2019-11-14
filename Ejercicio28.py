import os
import numpy as np
import matplotlib.pyplot as plt

plt.figure(1)

data = np.loadtxt("Ejercicio28.dat")

plt.plot(data[:,1], data[:,3])
plt.plot(data[:,5], data[:,6])

plt.title("Trajectoria de un proyectil con fricción y sin fricción")
plt.xlabel('X')
plt.ylabel('Y')
plt.savefig("Ejercicio28.png")
