from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv

AU = 1.5e11

t = np.array([]) #Create empty array for time
xsun = ysun = zsun = ssun = np.array([]) #Create arrays for Sun position in each dimension
xear = year = zear = sear = np.array([]) #Arrays for Earth
xjup = yjup = zjup = sjup = np.array([]) #Arrays for Jupiter
xmars = ymars = zmars = smars = np.array([]) #Arrays for Mars
xsat = ysat = zsat = ssat = np.array([]) #Arrays for Saturn
xura = yura = zura = sura = np.array([]) #Arrays for Uranus
xnep = ynep = znep = snep = np.array([]) #Arrays for Neptune

#sstar = t
#y = x**2
file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Sun_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsun = np.append(xsun,float(row[1]))
    ysun = np.append(ysun,float(row[2]))
    zsun = np.append(zsun,float(row[3]))
    ssun = np.append(ssun,float(row[4]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Earth_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    t = np.append(t,float(row[0]))
    xear = np.append(xear,float(row[1]))
    year = np.append(year,float(row[2]))
    zear = np.append(zear,float(row[3]))
    sear = np.append(sear,float(row[4]))
file.close()



file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Jupiter_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xjup = np.append(xjup,float(row[1]))
    yjup = np.append(yjup,float(row[2]))
    zjup = np.append(zjup,float(row[3]))
    sjup = np.append(sjup,float(row[4]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Mars_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xmars = np.append(xmars,float(row[1]))
    ymars = np.append(ymars,float(row[2]))
    zmars = np.append(zmars,float(row[3]))
    smars = np.append(smars,float(row[4]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Saturn_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsat = np.append(xsat,float(row[1]))
    ysat = np.append(ysat,float(row[2]))
    zsat = np.append(zsat,float(row[3]))
    ssat = np.append(ssat,float(row[4]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Uranus_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xura = np.append(xura,float(row[1]))
    yura = np.append(yura,float(row[2]))
    zura = np.append(zura,float(row[3]))
    sura = np.append(sura,float(row[4]))
file.close()

file = open('E:/Uni stuff not full of origin files/4th Year/Project/GitHub/PHY61008/Neptune_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xnep = np.append(xnep,float(row[1]))
    ynep = np.append(ynep,float(row[2]))
    znep = np.append(znep,float(row[3]))
    snep = np.append(snep,float(row[4]))
file.close()


fig = plt.figure(figsize =(20,20))
image = fig.add_subplot(111, projection = '3d')

#Set axis limits where all bodies should be visible
image.set_xlim(-40, 40)
image.set_ylim(-40, 40)
image.set_zlim(-40, 40)

plt.plot(xsun/AU, ysun/AU, zsun/AU, '-', c='yellow')
plt.plot(xear/AU, year/AU, zear/AU,'-', c='b')
plt.plot(xjup/AU, yjup/AU, zjup/AU, '-', c='saddlebrown')
plt.plot(xmars/AU, ymars/AU, zmars/AU,'-', c='r')
plt.plot(xsat/AU, ysat/AU, zsat/AU,'-', c='magenta')
plt.plot(xura/AU, yura/AU, zura/AU,'-', c='cyan')
plt.plot(xnep/AU, ynep/AU, znep/AU,'-', c='black')

plt.show()