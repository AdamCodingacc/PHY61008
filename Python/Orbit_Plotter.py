from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv

AU = 1.496e13 #in cm
pc = 3.086e18

t = np.array([]) #Create empty array for time
xsun = ysun = zsun = np.array([]) #Create arrays for Sun position in each dimension
xear = year = zear = np.array([]) #Arrays for Earth
xjup = yjup = zjup = np.array([]) #Arrays for Jupiter
xmars = ymars = zmars = np.array([]) #Arrays for Mars
xsat = ysat = zsat = np.array([]) #Arrays for Saturn
xura = yura = zura = np.array([]) #Arrays for Uranus
xnep = ynep = znep = np.array([]) #Arrays for Neptune


file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Sun_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Sun_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsun = np.append(xsun,float(row[1]))
    ysun = np.append(ysun,float(row[2]))
    zsun = np.append(zsun,float(row[3]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Earth_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Earth_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    t = np.append(t,float(row[0]))
    xear = np.append(xear,float(row[1]))
    year = np.append(year,float(row[2]))
    zear = np.append(zear,float(row[3]))
file.close()



file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Jupiter_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Jupiter_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xjup = np.append(xjup,float(row[1]))
    yjup = np.append(yjup,float(row[2]))
    zjup = np.append(zjup,float(row[3]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Mars_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Mars_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xmars = np.append(xmars,float(row[1]))
    ymars = np.append(ymars,float(row[2]))
    zmars = np.append(zmars,float(row[3]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Saturn_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Saturn_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsat = np.append(xsat,float(row[1]))
    ysat = np.append(ysat,float(row[2]))
    zsat = np.append(zsat,float(row[3]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Uranus_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Uranus_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xura = np.append(xura,float(row[1]))
    yura = np.append(yura,float(row[2]))
    zura = np.append(zura,float(row[3]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Neptune_Motion.csv')
#file = open('E:/Uni stuff not full of origin files/4th Year/Project/Neptune_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xnep = np.append(xnep,float(row[1]))
    ynep = np.append(ynep,float(row[2]))
    znep = np.append(znep,float(row[3]))
file.close()


fig = plt.figure(figsize =(20,20))
image = fig.add_subplot(111, projection = '3d')

"""
#Set axis limits where all bodies should be visible
image.set_xlim(-40, 40)
image.set_ylim(-40, 40)
image.set_zlim(-40, 40)
"""

plt.plot(xsun*pc/AU, ysun*pc/AU, zsun*pc/AU, '-', c='yellow')
plt.plot(xear*pc/AU, year*pc/AU, zear*pc/AU,'-', c='b')
plt.plot(xjup*pc/AU, yjup*pc/AU, zjup*pc/AU, '-', c='saddlebrown')
plt.plot(xmars*pc/AU, ymars*pc/AU, zmars*pc/AU,'-', c='r')
plt.plot(xsat*pc/AU, ysat*pc/AU, zsat*pc/AU,'-', c='magenta')
plt.plot(xura*pc/AU, yura*pc/AU, zura*pc/AU,'-', c='cyan')
plt.plot(xnep*pc/AU, ynep*pc/AU, znep*pc/AU,'-', c='black')

plt.show()