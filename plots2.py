import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d

#Variable initialization
x,y,z,error,dist,longitude,latitude=[],[],[],[],[],[],[]

#Put the data at the right place
var=''
reg_1=re.compile("-*.{12}[|]{1}")

with open("x_y_err_PEFRL.txt","r") as f:
    data=re.findall(reg_1,f.read())

for i in range(len(data)):
    for j in range(len(data[i])-1):
        var+=data[i][j]
    data[i]=float(var)
    var=''

for i in range(int(len(data)/7)):
    x+=[data[7*i]]    #Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y+=[data[7*i+1]]  #Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    z+=[data[7*i+2]]#Position on the z-absciss relative to the sun in the galactocentric frame in kpc
    error+=[data[7*i+3]]#Relative error on the energy of the pulsars in %
    dist+=[data[7*i+4]] #distance of the pulsars to the galactic center in kpc
    longitude+=[data[7*i+5]] #Longitude of the pulsars in degrees
    latitude+=[data[7*i+6]] #Latitude of the pulsars in degrees

#Positions plot
plt.figure(1)
plt.scatter(x,y,s=2,c='red',label='Simulation data')
plt.scatter([0],[0],c='yellow',marker='o',s=20,label='The Sun') #position of the sun
plt.xlim(-50,30)
plt.ylim(-20,25)
plt.title('Positions of the pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend()
plt.savefig('Positions_pulsars.png')

#Distance histogram
plt.figure(2)
plt.hist(dist,bins=201,range=(0,200),edgecolor='black',color='red',label='Simulation')
plt.legend()
plt.xlabel('d (kpc)')
plt.ylabel('Frequency')
plt.title('Distance from the galactic center')
plt.savefig('histo_dist_all_pulsars.png')

#Latitude histogram
plt.figure(3)
plt.hist(latitude,bins=121,range=(-60,60),edgecolor='black',color='red',label='Simulation')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
plt.title('Histogram of the latitude of all the pulsars')
plt.savefig('histo_latitude_all_pulsars.png')
