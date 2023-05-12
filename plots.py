import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha=[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out 
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
test_l=[]
Ba,B2=[],[]
B0=3.2e15 #constant to compute the surface magnetic field of the observations

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
    Ba+=[B0*(Pa[i]*P_dota[i])**0.5] #surface magnetic field in T

#sum_x,sum_y,sum_xy,sum_sqx,count=0,0,0,0,0
#for i in range(len(Pa)):
#    if P_dota[i]>0 and Ba[i]<5e8 and Ba[i]>5e6:
#        #if (Pa[i]> 3e-1 and P_dota[i]>1e-16):
#        #    if (Pa[i]<2 and P_dota[i]<1e-13):
#        sum_x+=np.log10(Pa[i])
#        sum_y+=np.log10(P_dota[i])
#        sum_xy+=np.log10(Pa[i])*np.log10(P_dota[i])
#        sum_sqx+=np.log10(Pa[i])**2
#        count+=1

#a_reglin=(sum_xy-((1/count))*(sum_x*sum_y))/(sum_sqx-(1/count)*(sum_x**2))
#b_reglin=(sum_y-a_reglin*sum_x)/count
#print(f'a : {a_reglin} \nb : {b_reglin} \n')

for i in range(len(P_dota)):
    if P_dota[i]!=0.0 and da[i]!=0.0 and E_dota[i]!=0:
        P2+=[Pa[i]]
        P_dot2+=[P_dota[i]]
        d2+=[da[i]]
        z2+=[za[i]]
        x2+=[xa[i]]
        y2+=[ya[i]]
        age2+=[agea[i]]
        E_dot2+=[E_dota[i]]
        B2+=[3.2e15*(Pa[i]*P_dota[i])**0.5]

with open("data_ATNF_for_c.txt","w") as f:
    for i in range(len(P2)):
        if B2[i]<5e8 and B2[i]>5e6:
            f.write(f"{P2[i]}")
            f.write("|")
            f.write(f"{P_dot2[i]}")
            f.write("|")
            f.write(f"{d2[i]}")
            f.write("|")
            f.write(f"{x2[i]}")
            f.write("|")
            f.write(f"{y2[i]}")
            f.write("|")
            f.write(f"{z2[i]}")
            f.write("|")
            f.write(f"{age2[i]}")
            f.write("|")
            f.write(f"{E_dot2[i]}")
            f.write("|")
            f.write("\n")

#count=0

#for i in range(len(z2)):
#    test=z2[i]/d2[i]
#    if test>=1 or test<=-1:
#        count+=1
#        print(i)
#    test_l+=[test]

#print(count)

for i in range(len(z2)):
    if d2[i]<1e-15 or i==9 or i==26:
        lat=0
        latitude2+=[lat]
    else:
        lat=np.arcsin(z2[i]/d2[i])*RAD
        latitude2+=[lat]
    r=(x2[i]**2+y2[i]**2)**0.5
    if x2[i]>=0:
        longi=np.arccos(-y2[i]/r)*RAD
        longitude+=[longi]
    else:
        longi=(np.arccos(-y2[i]/r)+np.pi)*RAD
        longitude+=[longi]

for i in range(len(age2)):
    log_age2+=[(np.log(age2[i]))/(np.log(10))]

for i in range(len(P2)):
    log_P2+=[(np.log(P2[i]))/(np.log(10))]
    log_Pdot2+=[(np.log(P_dot2[i]))/(np.log(10))]

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

for i in range(int(len(data)/11)):
    P+=[data[11*i]]    #Period of rotation of the pulsar in seconds
    P_dot+=[data[11*i+1]]  #Period derivative of rotation of the pulsar, no units 
    x+=[data[11*i+2]]#Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y+=[data[11*i+3]]#Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    age+=[data[11*i+4]] #Age of the pulsar in seconds 
    error+=[data[11*i+5]] #error on the positions of the pulsars (error computed with the energy)
    distance+=[data[11*i+6]] #Distance to the galactic center in kpc
    longitude+=[data[11*i+7]] #galactic latitude in degrees
    latitude+=[data[11*i+8]] #galactic longitude in degrees
    cos_alpha0+=[data[11*i+9]] #cosinus of the initial inclination angle
    cos_alpha+=[np.cos(data[11*i+10])] #cosinus of the inclination angle after all the evolution


for i in range(1,len(age)):
    logage=(np.log(age[i]/(365*24*60*60)))/(np.log(10))
    log_age+=[logage]

for i in range(len(P)):
    log_P+=[(np.log(P[i]))/(np.log(10))]
    log_Pdot+=[(np.log(P_dot[i]))/(np.log(10))]

#sum_x2,sum_y2,sum_xy2,sum_sqx2,count2=0,0,0,0,0
#for i in range(len(P)):
#    if P_dot[i]>0:
#        sum_x2+=np.log10(P[i])
#        sum_y2+=np.log10(P_dot[i])
#        sum_xy2+=np.log10(P[i])*np.log10(P_dot[i])
#        sum_sqx2+=np.log10(P[i])**2
#        count2+=1

#a_reglin2=(sum_xy2-((1/count2))*(sum_x2*sum_y2))/(sum_sqx2-(1/count2)*(sum_x2**2))
#b_reglin2=(sum_y2-a_reglin2*sum_x2)/count2
#print(f'a_simu : {a_reglin2} \n b_simu : {b_reglin2} \n')

#Lists depending on the pulsar emission type
P_radio,P_dot_radio,x_radio,y_radio,age_radio,error_radio,distance_radio=[],[],[],[],[],[],[]
P_gamma,P_dot_gamma,x_gamma,y_gamma,age_gamma,error_gamma,distance_gamma=[],[],[],[],[],[],[]
P_radio_gamma,P_dot_radio_gamma,x_radio_gamma,y_radio_gamma,age_radio_gamma,error_radio_gamma,distance_radio_gamma=[],[],[],[],[],[],[]
for i in range(len(P)):
    if type_pulsar[i]==1:
        P_radio+=[P[i]]
        P_dot_radio+=[P_dot[i]]
        x_radio+=[x[i]]
        y_radio+=[y[i]]
        age_radio+=[age[i]]
        distance_radio+=[distance[i]]
    elif type_pulsar[i]==2:
        P_gamma+=[P[i]]
        P_dot_gamma+=[P_dot[i]]
        x_gamma+=[x[i]]
        y_gamma+=[y[i]]
        age_gamma+=[age[i]]
        distance_gamma+=[distance[i]]
    elif type_pulsar[i]==3:
        P_radio_gamma+=[P[i]]
        P_dot_radio_gamma+=[P_dot[i]]
        x_radio_gamma+=[x[i]]
        y_radio_gamma+=[y[i]]
        age_radio_gamma+=[age[i]]
        distance_radio_gamma+=[distance[i]]

#Plot the death line of the article from Mitra et al. (2019)
T_6=2
eta=0.15
alpha_l=45*np.pi/180
b=40
const=(3.16e-4*T_6*1e-15)/((eta)**2*b*(np.cos(alpha_l))**2)
P_line=[i for i in np.arange(1e-2,1e1,0.001)]
Pdot_line=[const*(i**2) for i in np.arange(1e-2,1e1,0.001)]

#L1=[i for i in np.arange(1e-2,1e1,0.1)]
#L2=[10**(a_reglin*np.log10(i)+b_reglin) for i in np.arange(1e-2,1e1,0.1)]
#L3=[10**(a_reglin2*np.log10(i)+b_reglin2) for i in np.arange(1e-2,1e1,0.1)]
#Make the plots
#P-Pdot plot all pulsars
plt.figure(1)
plt.scatter(P,P_dot,c='red',marker='o',s=5,label='Simulation data')
plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data')
plt.plot(P_line,Pdot_line)
#plt.plot(L1,L2)
#plt.plot(L1,L3)
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
plt.scatter(x,y,s=2,c='red',label='Simulation data')
plt.scatter(x2,y2,s=2,c='blue',label='ATNF data')
plt.scatter([0],[8.5],c='yellow',marker='o',s=20,label='The Sun') #position of the sun
plt.xlim(-50,30)
plt.ylim(-30,35)
plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend()
plt.savefig('Positions_detected_pulsars.png')

#Distance histogram
plt.figure(6)
plt.hist(distance,bins=70,range=(0,25),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.hist(d2,bins=70,range=(0,25),edgecolor='black',color='blue',alpha=0.5,label='ATNF data')
plt.legend()
plt.xlabel('d (kpc)')
plt.ylabel('Frequency')
plt.title('Distance to earth of the pulsars')
plt.savefig('histo_dist.png')

#Age histogram
plt.figure(7)
plt.hist(log_age,bins=130,range=(1,11),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.hist(log_age2,bins=130,range=(1,11),edgecolor='black',color='blue',alpha=0.5,label='ATNF data')
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
plt.title('Histogram of the age of the detected pulsars')
plt.savefig('histo_age.png')

#Latitude histogram
plt.figure(8)
plt.hist(latitude,bins=121,range=(-60,60),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.hist(latitude2,bins=121,range=(-60,60),edgecolor='black',color='blue',alpha=0.5,label='ATNF data')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
plt.title('Histogram of the latitude of the detected pulsars')
plt.savefig('histo_latitude.png')

#Log(P) histogram
plt.figure(9)
plt.hist(log_P,bins=150,range=(-2,1.5),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.hist(log_P2,bins=150,range=(-2,1.5),edgecolor='black',color='blue',alpha=0.5,label='ATNF data')
plt.legend()
plt.xlabel('Log(P) (P in s)')
plt.ylabel('Frequency')
plt.title('Histogram of the rotation period of the detected pulsars')
plt.savefig('histo_period.png')

#Log(P_dot) histogram
plt.figure(10)
plt.hist(log_Pdot,bins=150,range=(-20,-10),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.hist(log_Pdot2,bins=150,range=(-20,-10),edgecolor='black',color='blue',alpha=0.5,label='ATNF data')
plt.legend()
plt.xlabel('Log(Pdot)')
plt.ylabel('Frequency')
plt.title('Histogram of the rotation period derivative of the detected pulsars')
plt.savefig('histo_pdot.png')

#cos(alpha0) and cos(alpha) histogram
plt.figure(11)
plt.hist(cos_alpha,bins=100,range=(0,1),edgecolor='black',color='red',alpha=0.5,label='cos(alpha)')
plt.hist(cos_alpha0,bins=100,range=(0,1),edgecolor='black',color='blue',alpha=0.5,label='cos(alpha0)')
plt.legend()
plt.xlabel('cos(alpha) and cos(alpha0)')
plt.ylabel('Frequency')
plt.title('Histogram of the inclination angle of the detected pulsars')
plt.savefig('histo_cosalpha.png')
