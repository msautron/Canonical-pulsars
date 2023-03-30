import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar=[],[],[],[],[],[],[]
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[]
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2=[],[],[],[],[],[],[],[]

#Put the data at the right place
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("[-+]?\d*[.]\d*[Ee]*[-+]*\d*")

with open("ATNF_data.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data=re.findall(reg_1,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data_type=re.findall(reg_2,f.read())

#Get the data of the ATNF catalogue
for i in range(int(len(data2)/8)):
    Pa+=[float(data2[i*8])] #Period of the rotation of the pulsar in seconds
    P_dota+=[float(data2[i*8+1])] #Period derivative of the rotation of the pulsar no units
    da+=[float(data2[i*8+2])] #Distance to us in kpc
    za+=[float(data2[i*8+3])] #Z position in the galactocentric frame in kpc
    xa+=[float(data2[i*8+4])] #X position in the galactocentric frame in kpc
    ya+=[float(data2[i*8+5])] #age of the pulsar in years
    agea+=[float(data2[i*8+6])] #age of the pulsars in seconds
    E_dota+=[float(data2[i*8+7])]  #Spin down power of the pulsar in ergs/s

for i in range(len(P_dota)):
    if P_dota[i]!=0.0:
        P2+=[Pa[i]]
        P_dot2+=[P_dota[i]]
        d2+=[da[i]]
        z2+=[0.015-za[i]]
        x2+=[xa[i]]
        y2+=[8.5-ya[i]]
        age2+=[agea[i]]
        E_dot2+=[E_dota[i]]

#get the data of the P_Pdot_positions file (simulation)
for i in range(len(data_type)):
    var2+=data_type[i][1]
    type_pulsar+=[int(var2)]
    var2=''


for i in range(len(data)):
    for j in range(len(data[i])-1):
        var+=data[i][j]
    data[i]=float(var)
    var=''

for i in range(int(len(data)/6)):
    P+=[data[6*i]]    #Period of rotation of the pulsar in seconds
    P_dot+=[data[6*i+1]]  #Period derivative of rotation of the pulsar, no units 
    x+=[data[6*i+2]]#Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y+=[data[6*i+3]]#Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    age+=[data[6*i+4]] #Age of the pulsar in seconds 
    error+=[data[6*i+5]] #error on the positions of the pulsars (error computed with the energy)

#Lists depending on the pulsar emission type
P_radio,P_dot_radio,x_radio,y_radio,age_radio,error_radio=[],[],[],[],[],[]
P_gamma,P_dot_gamma,x_gamma,y_gamma,age_gamma,error_gamma=[],[],[],[],[],[]
P_radio_gamma,P_dot_radio_gamma,x_radio_gamma,y_radio_gamma,age_radio_gamma,error_radio_gamma=[],[],[],[],[],[]
for i in range(len(P)):
    if type_pulsar[i]==1:
        P_radio+=[P[i]]
        P_dot_radio+=[P_dot[i]]
        x_radio+=[x[i]]
        y_radio+=[y[i]]
        age_radio+=[age[i]]
    elif type_pulsar[i]==2:
        P_gamma+=[P[i]]
        P_dot_gamma+=[P_dot[i]]
        x_gamma+=[x[i]]
        y_gamma+=[y[i]]
        age_gamma+=[age[i]]
    elif type_pulsar[i]==3:
        P_radio_gamma+=[P[i]]
        P_dot_radio_gamma+=[P_dot[i]]
        x_radio_gamma+=[x[i]]
        y_radio_gamma+=[y[i]]
        age_radio_gamma+=[age[i]]

#Make the plots
#P-Pdot plot all pulsars
plt.figure(1)
plt.scatter(P,P_dot,c='red',marker='o',s=10,label='Simulation data')
plt.scatter(P2,P_dot2,c='blue',marker='o',s=10,label='ATNF data')
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
plt.title("Spin period derivative - Spin period diagram")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.legend()
plt.savefig('P_Pdot_plot.png')

#P-Pdot plot radio pulsars only 
plt.figure(2)
plt.scatter(P_radio,P_dot_radio,c='red',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
plt.title("Spin period derivative - Spin period diagram for radio pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_r.png')

#P-Pdot plot gamma pulsars only
plt.figure(3)
plt.scatter(P_gamma,P_dot_gamma,c='green',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
plt.title("Spin period derivative - Spin period diagram for gamma pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_g.png')

#P-Pdot plot radio gamma pulsars
plt.figure(4)
plt.scatter(P_radio_gamma,P_dot_radio_gamma,c='purple',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
plt.title("Spin period derivative - Spin period diagram for radio gamma pulsars")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_rg.png')

#Positions plot
plt.figure(5)
plt.scatter(x,y,s=10,c='red',label='Simulation data')
plt.scatter(x2,y2,s=10,c='blue',label='ATNF data')
plt.scatter([0],[0],c='yellow',marker='o',s=40) #position of the sun
plt.xlim(-50,30)
plt.ylim(-20,25)
plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.savefig('Positions_detected_pulsars.png')

#Tests
