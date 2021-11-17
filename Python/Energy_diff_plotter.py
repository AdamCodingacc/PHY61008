from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv

E_diff = E_percent = t = np.array([])

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Energy_err.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    t = np.append(t,float(row[0]))
    E_diff = np.append(E_diff,float(row[1]))
    E_percent = np.append(E_percent,float(row[2]))
file.close()

#plt.plot(E_diff, t)
plt.plot(E_percent, t)