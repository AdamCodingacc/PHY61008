from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import csv

AU = 1.5e11

t = np.array([]) #Create empty array for time
xsun, ysun, zsun, ssun = t #Create arrays for Sun position in each dimension
xear, year, zear, sear = t #Arrays for Earth
xjup, yjup, zjup, sjup = t #Arrays for Jupiter
xmars, ymars, zmars, smars = t #Arrays for Mars
xsat, ysat, zsat, ssat = t #Arrays for Saturn
xura, yura, zura, sura = t #Arrays for Uranus
xnep, ynep, znep, snep = t #Arrays for Neptune

#sstar = t
#y = x**2
file = open('F:/Fortran/Projects/3 Body Problem/Sun_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsun = np.append(xsun,float(row[1]))
    ysun = np.append(ysun,float(row[2]))
    zsun = np.append(zsun,float(row[3]))
    ssun = np.append(ssun,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Earth_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    t = np.append(t,float(row[0]))
    xear = np.append(xear,float(row[1]))
    year = np.append(year,float(row[2]))
    zear = np.append(zear,float(row[3]))
    sear = np.append(sear,float(row[4]))
file.close()



file = open('F:/Fortran/Projects/3 Body Problem/Jupiter_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xjup = np.append(xjup,float(row[1]))
    yjup = np.append(yjup,float(row[2]))
    zjup = np.append(zjup,float(row[3]))
    sjup = np.append(sjup,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Mars_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xmars = np.append(xmars,float(row[1]))
    ymars = np.append(ymars,float(row[2]))
    zmars = np.append(zmars,float(row[3]))
    smars = np.append(smars,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Saturn_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xsat = np.append(xsat,float(row[1]))
    ysat = np.append(ysat,float(row[2]))
    zsat = np.append(zsat,float(row[3]))
    ssat = np.append(ssat,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Uranus_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xura = np.append(xura,float(row[1]))
    yura = np.append(yura,float(row[2]))
    zura = np.append(zura,float(row[3]))
    sura = np.append(sura,float(row[4]))
file.close()

file = open('F:/Fortran/Projects/3 Body Problem/Neptune_Motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    xnep = np.append(xnep,float(row[1]))
    ynep = np.append(ynep,float(row[2]))
    znep = np.append(znep,float(row[3]))
    snep = np.append(snep,float(row[4]))
file.close()

"""
file = open('F:/Fortran/Projects/3 Body Problem/star motion.csv')
look = csv.reader(file,delimiter = ',')
for row in look:
    sstar = np.append(sstar,float(row[1]))
file.close()
"""






fig = plt.figure(figsize =(20,20))
axis = fig.add_subplot(111, projection = '3d',)
#plt.axis("off")
limit= 20
axis.set_xlim(-1*limit,limit)
axis.set_ylim(-1*limit,limit)
axis.set_zlim(-1*limit,limit)

plt.plot(xsun/AU, ysun/AU, zsun/AU, c='yellow')
plt.plot(xear/AU,year/AU,zear/AU,'-',c='b')
plt.plot(xjup/AU, yjup/AU, zjup/AU, '-',c='orange')
plt.plot(xmars/AU,ymars/AU,zmars/AU,'-',c='r')
plt.plot(xsat/AU,ysat/AU,zsat/AU,'-',c='magenta')
plt.plot(xura/AU,yura/AU,zura/AU,'-',c='cyan')
plt.plot(xnep/AU,ynep/AU,znep/AU,'-',c='purple')
#plt.plot(2*x/AU,y/AU,z/AU, c = 'r')
#axis.view_init(0,-90)
plt.show()
#plt.plot(t,se/AU)
#print(t)
#print(x)
#print(row)

"""
print(max(se/AU),min(se/AU))
print(max(sj/AU),min(sj/AU))
print(max(sm/AU),min(sm/AU))
print(max(sst/AU),min(sst/AU))
print(max(su/AU),min(su/AU))
print(max(sn/AU),min(sn/AU))
ee = (max(se/AU) - min(se/AU))/(max(se/AU) + min(se/AU))
ej = (max(sj/AU) - min(sj/AU))/(max(sj/AU) + min(sj/AU))
em = (max(sm/AU) - min(sm/AU))/(max(sm/AU) + min(sm/AU))
est = (max(sst/AU) - min(sst/AU))/(max(sst/AU) + min(sst/AU))
eu = (max(su/AU) - min(su/AU))/(max(su/AU) + min(su/AU))
en = (max(sn/AU) - min(sn/AU))/(max(sn/AU) + min(sn/AU))
Eccentricities = [ee,ej,em,est,eu,en]


i = 1
for e in Eccentricities:
    if e > 0.1:
        print("Perturbation Detected")
        print(i,e)
    i += 1
i = 1
print(min(sstar)/AU)

"""
#%%
#fig,axis = plt.subplots(figsize = (20,20))
#plt.plot(t/(365*24),se/AU)
#axis.set_ylim(,0.9974)
#plt.plot(t/365,sj/AU)
#plt.show()