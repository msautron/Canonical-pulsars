import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out 
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
test_l=[]
Ba,B2=[],[]
B0=3.2e15 #constant to compute the surface magnetic field of the observations
P_selected,Pdot_selected=[],[]
#B_init=[]

#Put the data at the right place
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("[-+]?\d*[.]\d*[Ee]*[-+]*\d*")
reg_4=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
#reg_5=re.compile("-*.{12}")

with open("ATNF_data.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data=re.findall(reg_1,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data_type=re.findall(reg_2,f.read())

with open("xy_coord_dirson22.dat","r") as f :
    data1_dirson22=re.findall(reg_4,f.read())

with open("glat_dirson22.dat","r") as f :
    data2_dirson22=re.findall(reg_4,f.read())

with open("dist_dirson22.dat","r") as f:
    data3_dirson22=re.findall(reg_4,f.read())

#Get the data from Dirson et al. (2022) for xy coord and latitude
x_dirson22,y_dirson22,lat_dirson22,dist_dirson22=[],[],[],[]
for i in range(int(len(data1_dirson22)/2)):
    x_dirson22.append(float(data1_dirson22[2*i]))
    y_dirson22.append(float(data1_dirson22[2*i+1]))

for i in range(int(len(data2_dirson22))):
    lat_dirson22.append(float(data2_dirson22[i]))

for i in range(int(len(data3_dirson22))):
    dist_dirson22.append(float(data3_dirson22[i]))

#with open("verif_lognorm.txt","r") as f:
#    data_lognorm=re.findall(reg_4,f.read())

#Get the data to verify it is a log normal distribution that we are using
#for i in range(len(data_lognorm)):
#    B_init+=[float(data_lognorm[i])]

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

for i in range(len(P_dota)):
    #P2.append(Pa[i])
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

d_selec,z_selec,x_selec,y_selec,age_selec,E_dot_selec,B_selec=[],[],[],[],[],[],[]
d_magne,z_magne,x_magne,y_magne,age_magne,E_dot_magne,B_magne,P_magne,Pdot_magne=[],[],[],[],[],[],[],[],[]
P_ms,Pdot_ms,d_ms,z_ms,x_ms,y_ms,age_ms,E_dot_ms,B_ms=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dot2)): #Selection of P and P_dot to get only the pulsars for the optimization program
    if B2[i]<4.4e9 and B2[i]>6e6:
        P_selected+=[P2[i]]
        Pdot_selected+=[P_dot2[i]]
        d_selec.append(d2[i])
        z_selec.append(z2[i])
        x_selec.append(x2[i])
        y_selec.append(y2[i])
        age_selec.append(age2[i])
        E_dot_selec.append(E_dot2[i])
        B_selec.append(B2[i])
    elif B2[i]>=4.4e9:
        P_magne+=[P2[i]]
        Pdot_magne+=[P_dot2[i]]
        d_magne.append(d2[i])
        z_magne.append(z2[i])
        x_magne.append(x2[i])
        y_magne.append(y2[i])
        age_magne.append(age2[i])
        E_dot_magne.append(E_dot2[i])
        B_magne.append(B2[i])
    elif B2[i]<=6e6:
        P_ms+=[P2[i]]
        Pdot_ms+=[P_dot2[i]]
        d_ms.append(d2[i])
        z_ms.append(z2[i])
        x_ms.append(x2[i])
        y_ms.append(y2[i])
        age_ms.append(age2[i])
        E_dot_ms.append(E_dot2[i])
        B_ms.append(B2[i])

with open("data_ATNF_for_c.txt","w") as f:
    for i in range(len(P2)):
        if B2[i]<4.4e9 and B2[i]>6e6:
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

#Computation of latitude and longitude for the ATNF data
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

latitude_selec,longitude_selec=[],[]
for i in range(len(z_selec)):
    if d_selec[i]<1e-15:
        lat=0
        latitude_selec+=[lat]
    else:
        lat=np.arcsin(z_selec[i]/d_selec[i])*RAD
        latitude_selec+=[lat]
    r=(x_selec[i]**2+y_selec[i]**2)**0.5
    if x_selec[i]>=0:
        longi=np.arccos(-y_selec[i]/r)*RAD
        longitude_selec+=[longi]
    else:
        longi=(np.arccos(-y_selec[i]/r)+np.pi)*RAD
        longitude_selec+=[longi]
#Log calculation for age, P and Pdot for the ATNF data
log_age_selec=[]
for i in range(len(age2)):
    log_age2+=[(np.log(age2[i]))/(np.log(10))]

for i in range(len(age_selec)):
    log_age_selec+=[np.log10(age_selec[i])]

log_P_selec,log_Pdot_selec=[],[]
for i in range(len(P_selected)):
    log_P_selec+=[np.log10(P_selected[i])]
    log_Pdot_selec+=[np.log10(Pdot_selected[i])]

log_P_ms,log_Pdot_ms=[],[]
for i in range(len(P_ms)):
    log_P_ms+=[np.log10(P_ms[i])]
    log_Pdot_ms+=[np.log10(Pdot_ms[i])]

log_P_magne,log_Pdot_magne=[],[]
for i in range(len(P_magne)):
    log_P_magne+=[np.log10(P_magne[i])]
    log_Pdot_magne+=[np.log10(Pdot_magne[i])]

log_Pa=[]
for i in range(len(Pa)):
    log_Pa+=[np.log10(Pa[i])]

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

for i in range(int(len(data)/20)):
    P+=[data[20*i]]    #Period of rotation of the pulsar in seconds
    P_dot+=[data[20*i+1]]  #Period derivative of rotation of the pulsar, no units 
    x+=[data[20*i+2]]#Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y+=[data[20*i+3]]#Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    age+=[data[20*i+4]] #Age of the pulsar in seconds 
    error+=[data[20*i+5]] #error on the positions of the pulsars (error computed with the energy)
    distance+=[data[20*i+6]] #Distance to the galactic center in kpc
    longitude+=[data[20*i+7]] #galactic latitude in degrees
    latitude+=[data[20*i+8]] #galactic longitude in degrees
    cos_alpha0+=[data[20*i+9]] #cosinus of the initial inclination angle
    cos_alpha+=[np.cos(data[20*i+10])] #cosinus of the inclination angle after all the evolution
    Bf+=[data[20*i+11]] #Magnetic field of the pulsar today, in tesla
    z+=[data[20*i+12]] #Position on the z-absciss relative to the sun in the galactocentric frame in kpc
    vx+=[data[20*i+13]] #Velocity on the x-absciss of the pulsar today, in km/s
    vy+=[data[20*i+14]] #Velocity on the y-absciss of the pulsar today, in km/s
    vz+=[data[20*i+15]] #Velocity on the z-absciss of the pulsar today, in km/s
    vx0+=[data[20*i+16]] #Velocity on the x-absciss of the pulsar initially, in km/s
    vy0+=[data[20*i+17]] #Velocity on the y-absciss of the pulsar initially, in km/s
    vz0+=[data[20*i+18]] #Velocity on the z-absciss of the pulsar initially, in km/s
    PA+=[(180/np.pi)*np.arccos(data[20*i+19])]

#Edot computation 
Edot=[]
Inertia=1e38
count_rad,count_gam,count_radgam,count_radgam_E_big,count_gam_E_big,count_rad_E_big,count_radgam_E_bigbig,count_gam_E_bigbig,count_rad_E_bigbig=0,0,0,0,0,0,0,0,0
for i in range(len(P)):
    Edot+=[4*np.pi**2*Inertia*P_dot[i]*(P[i]**(-3))]

for i in range(0,len(age)):
    logage=(np.log(age[i]/(365*24*60*60)))/(np.log(10))
    log_age+=[logage]

for i in range(len(P)):
    log_P+=[(np.log(P[i]))/(np.log(10))]
    log_Pdot+=[(np.log(P_dot[i]))/(np.log(10))]

#Prep death line Ruderman & Sutherland 1975
R_NS=12000
mu_0=1.25663706212e-6 
c_light=2.997924858e8
P_dot_death,P_dot_death2=[],[]
P_death=[np.log10(i) for i in np.arange(1e-2,3e1,0.01)]
P_death2=[10**(P_death[i]) for i in range(len(P_death))]
for i in range(len(P_death)):
    P_dot_death+=[3*P_death[i]+np.log10(((16*(np.pi**3)*(R_NS**6)*(1+((np.sin(45*np.pi/180)**2))))*(17000000**2))/(Inertia*mu_0*(c_light**3)))]
    P_dot_death2+=[10**(P_dot_death[i])]

#B(P,Pdot) computation in order to compare with the decaying Bf
B_ppdot,v_norm,v0_norm,err_rel_B=[],[],[],[]
for i in range(len(P)):
    B_ppdot+=[((Inertia*mu_0*(c_light**3)*P_dot[i]*P[i])/(16*(np.pi**3)*(R_NS**6)*(1+(np.sin(np.arccos(cos_alpha[i]))**2))))**0.5]
    v_norm+=[(vx[i]**2+vy[i]**2+vz[i]**2)**0.5]
    v0_norm+=[(vx0[i]**2+vy0[i]**2+vz0[i]**2)**0.5]
    err_rel_B+=[(np.abs(B_ppdot[i]-Bf[i]))/np.abs(B_ppdot[i])]

#Prep plot period old pulsars and plot old pulsars spin-velocity angle
P_old=[]
PA_old,PA_young=[],[]
for i in range(len(log_age)):
    if log_age[i]>7.5 and log_age[i]<9:
        P_old+=[log_P[i]]

for i in range(len(PA)):
    if log_age[i]>=7:
        PA_old.append(PA[i])
    elif log_age[i]<7:
        PA_young.append(PA[i])


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

#d3=[]
#for i in range(len(d2)):
#    if B2[i]>=1e6 and B2[i]<=1e9:
#        d3.append(d2[i])

#print(len(d3))
#print(len(d2))

#Plot the death line of the article from Mitra et al. (2019)
T_6=2
T_6_max=2.8
T_6_min=1.9
eta=0.15
alpha_l=45*np.pi/180
alpha_l_max=65*np.pi/180
alpha_l_min=0*np.pi/180
b=40
b_min=60
b_max=30
const=(3.16e-4*(T_6**4)*1e-15)/((eta)**2*b*(np.cos(alpha_l))**2)
const_min=(3.16e-4*(T_6_min**4)*1e-15)/((eta)**2*b_min*(np.cos(alpha_l_min))**2)
const_max=(3.16e-4*(T_6_max**4)*1e-15)/((eta)**2*b_max*(np.cos(alpha_l_max))**2)
P_line=[10**(np.log10(i)) for i in np.arange(1e-3,3e1,0.01)]
Pdot_line=[10**np.log10(const*(i**2)) for i in np.arange(1e-3,3e1,0.01)]
Pdot_line4=[Pdot_line[i]*10**(-0.55) for i in range(len(Pdot_line))]
Pdot_line5=[Pdot_line[i]*10**(1.15) for i in range(len(Pdot_line))]
Pdot_line2=[10**np.log10(const_min*(i**2)) for i in np.arange(1e-3,3e1,0.01)]
Pdot_line3=[10**np.log10(const_max*(i**2)) for i in np.arange(1e-3,3e1,0.01)]

#Plot the death line of the article of Chen & Ruderman (1993)
const_CR93=10**((43.8-16)/2)*16*(np.pi**3)*(R_NS**6)*(1+(np.sin(alpha_l)**2))/(Inertia*mu_0*(c_light**3))
P_dot_line_CR93=[10**np.log10(const_CR93*(i**(2))) for i in np.arange(1e-3,3e1,0.01)]

#Plot the line with the critical magnetic field 4.4e9 T, distinguishing magnetar from canonical pulsars
B_crit=4.4e9
B_linecrit=[10**np.log10(((B_crit/B0)**2)*(i**(-1))) for i in np.arange(1e-3,3e1,0.01)]


#Check if the pulsars are really acceptable with the ratio B/P^2 - 0.17e8 or with the Pdot of the death of Mitra et al. (2019)
Diff=[]
count_abno=0
for i in range(len(P)):
    #Diff+=[(Bf[i]/((P[i])**2))]
    #Diff+=[(Bf[i]**2)/((P[i])**3)]
    #Diff+=[np.log10(Bf[i]) - 2*np.log10(P[i]) - np.log10(0.17) -8.0]
    Diff+=[(3.16e-4*(T_6**4)*1e-15*(P[i]**2))/((eta)**2*b*(np.cos(alpha_l))**2)]
    if Diff[i]>P_dot[i]:
    #if Diff[i]>10**(13.9):
        count_abno+=1
        #print(i)
        #print(Diff[i])
        #print(P[i])
print(f"The number of pulsars which should be dead is : {count_abno}\nlength list :{len(Diff)}\n")
#print(Diff[0])
#print(Bf[0])
#print(P[0])

#Remove all the points below the death line 
#P_dot_line=[]
#nb_rm=0
#for i in range(len(P)):
    #P_dot_line+=[3*np.log10(P[i])+np.log10(((16*(np.pi**3)*(R_NS**6)*(1+((np.sin(45*np.pi/180)**2))))*(17000000**2))/(Inertia*mu_0*(c_light**3)))]
    #P_dot_line+=[(3.16e-4*(T_6**4)*1e-15*(P[i]**2))/((eta)**2*b*(np.cos(alpha_l))**2)]

#for j in range(len(P_dot_line)):
    #if P_dot_line[j]>P_dot[j-nb_rm]:
        #del P[j-nb_rm]
        #del P_dot[j-nb_rm]
        #del x[j-nb_rm] #Position on the x-absciss relative to the sun in the galactocentric frame in kpc
        #del y[j-nb_rm] #Position on the y-absciss relative to the sun in the galactocentric frame in kpc
        #del age[j-nb_rm] #Age of the pulsar in seconds
        #del error[j-nb_rm] #error on the positions of the pulsars (error computed with the energy)
        #del distance[j-nb_rm] #Distance to the galactic center in kpc
        #del longitude[j-nb_rm] #galactic latitude in degrees
        #del latitude[j-nb_rm] #galactic longitude in degrees
        #del cos_alpha0[j-nb_rm] #cosinus of the initial inclination angle
        #del cos_alpha[j-nb_rm] #cosinus of the inclination angle after all the evolution
        #del Bf[j-nb_rm] #Magnetic field of the pulsar today, in tesla
        #del z[j-nb_rm] #Position on the z-absciss relative to the sun in the galactocentric frame in kpc
        #del vx[j-nb_rm] #Velocity on the x-absciss of the pulsar today, in km/s
        #del vy[j-nb_rm] #Velocity on the y-absciss of the pulsar today, in km/s
        #del vz[j-nb_rm] #Velocity on the z-absciss of the pulsar today, in km/s
        #del vx0[j-nb_rm] #Velocity on the x-absciss of the pulsar initially, in km/s
        #del vy0[j-nb_rm] #Velocity on the y-absciss of the pulsar initially, in km/s
        #del vz0[j-nb_rm]
        #del Edot[j-nb_rm]
        #del type_pulsar[j-nb_rm]
        #nb_rm+=1

#Display of the statistics on all pulsar
for i in range(len(Edot)):
    if Edot[i] > 1e31 and type_pulsar[i]==1:
        count_rad_E_bigbig+=1
    if Edot[i] > 1e31 and type_pulsar[i]==2:
        count_gam_E_bigbig+=1
    if Edot[i] > 1e31 and type_pulsar[i]==3:
        count_radgam_E_bigbig+=1
    if Edot[i] > 1e28 and type_pulsar[i]==1:
        count_rad_E_big+=1
    if Edot[i] > 1e28 and type_pulsar[i]==2:
        count_gam_E_big+=1
    if Edot[i] > 1e28 and type_pulsar[i]==3:
        count_radgam_E_big+=1
    if type_pulsar[i]==1:
        count_rad+=1
    if type_pulsar[i]==2:
        count_gam+=1
    if type_pulsar[i]==3:
        count_radgam+=1

print(f"Number of radio pulsars with Edot > 1e31 W : {count_rad_E_bigbig}\nNumber of gamma pulsars with Edot > 1e31 W : {count_gam_E_bigbig}")
print(f"Number of radio-gamma pulsars with Edot > 1e31 W : {count_radgam_E_bigbig}")
print(f"Number of radio pulsars with Edot > 1e28 W : {count_rad_E_big}")
print(f"Number of gamma pulsars with Edot > 1e28 W : {count_gam_E_big}")
print(f"Number of radio-gamma pulsars with Edot > 1e28 W : {count_radgam_E_big}")
print(f"Number of radio pulsars : {count_rad}")
print(f"Number of gamma pulsars : {count_gam}")
print(f"Number of radio-gamma pulsars : {count_radgam}")

histSIM,xsim_edges,ysim_edges=np.histogram2d(log_P,log_Pdot,bins=(30,30))
histobs,xobs_edges,yobs_edges=np.histogram2d(log_P_selec,log_Pdot_selec,bins=(30,30))
histSIM=np.rot90(histSIM)
#histSIM=np.rot90(histSIM)
histobs=np.rot90(histobs)
#histobs=np.rot90(histobs)

hist_diff=((histSIM/len(log_P))-(histobs/len(log_P_selec)))/(((histSIM/len(log_P))+(histobs/len(log_P_selec)))**1.0)
chi2_tab=(histSIM-histobs)**2/histobs
chi2_tab=np.where(np.isfinite(chi2_tab), chi2_tab, 0)
chi2=np.nansum(chi2_tab)
print(chi2)
print(np.shape(chi2_tab))

#Make the plots
#test
#plt.figure(12)
#plt.scatter(P_selected,Pdot_selected,c='blue',marker='o',s=5,label='ATNF data')
#plt.scatter(log_P,log_Pdot,c='red',marker='o',s=5,label='Simulation data')
#plt.xlim(1e-2,1e1)
#plt.ylim(1e-20,1e-10)
#plt.yscale('log')
#plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
#plt.xlabel('Spin period s')
#plt.ylabel('Spin period derivative s.s^-1')
#plt.legend()
#plt.savefig('selected_P_Pdot.png')

#print(len(P_dot2))
#print(len(Pa))

condition = [Pdot2 < Pdot3 for Pdot2, Pdot3 in zip(Pdot_line2,Pdot_line3)]

#P-Pdot plot only canonical population
plt.figure(30)
plt.scatter(P_selected,Pdot_selected,c='blue',marker='o',s=10,label='ATNF data')
plt.legend()
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_selected.png')
plt.close()

#P-Pdot plot all pulsars
plt.figure(1)
plt.scatter(P,P_dot,c='red',marker='o',s=5,label='Simulation data',zorder=2)
#plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data') #whole pop
plt.scatter(P_selected,Pdot_selected,c='blue',marker='o',s=5,label='ATNF data : canonical pulsars') #Only canonical pop
plt.scatter(P_magne,Pdot_magne,c='orange',marker='s',s=5,label='ATNF data : magnetars') #Only magnetars pop
plt.scatter(P_ms,Pdot_ms,c='purple',marker='^',s=5,label='ATNF data : millisecond pulsars') #Only ms pop
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
plt.plot(P_line,Pdot_line3,c='brown')
plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )
plt.xlim(1e-2,3e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel('spin period (s)')
plt.ylabel('derivative of the spin period (s.s^-1)')
plt.legend(fontsize='x-small')
plt.savefig('P_Pdot_plot.png',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(31)
#plt.figure(figsize=(6,6))
#plt.imshow((histSIM-histobs)/((histSIM+histobs)**0.5),extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(xsim_edges[:-1],ysim_edges[:-1],bins=(xsim_edges,ysim_edges),cmap='RdBu',weights=hist_diff.flatten())
plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel('Log P (P in s)')
plt.ylabel('Log Pdot (Pdot in s.s^-1)')
plt.savefig('2D P_Pdot_plot_sim.png',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(32)
#plt.figure(figsize=(6,6))
#plt.imshow((histSIM-histobs)/((histSIM+histobs)**0.5),extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.hist2d(log_P_selec,log_Pdot_selec,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel('Log P (P in s)')
plt.ylabel('Log Pdot (Pdot in s.s^-1)')
plt.savefig('2D P_Pdot_plot_obs.png',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(33)
#plt.figure(figsize=(6,6))
plt.imshow(hist_diff,extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(xsim_edges[:-1],ysim_edges[:-1],bins=(xsim_edges,ysim_edges),cmap='RdBu',weights=hist_diff.flatten())
#plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel('Log P (P in s)')
plt.ylabel('Log Pdot (Pdot in s.s^-1)')
plt.savefig('2D P_Pdot_plot.png',dpi=300)
plt.close()

#P-Pdot plot radio pulsars only 
plt.figure(2)
plt.scatter(P_radio,P_dot_radio,c='red',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_r.png')
plt.close()

#P-Pdot plot gamma pulsars only
plt.figure(3)
plt.scatter(P_gamma,P_dot_gamma,c='green',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for gamma pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_g.png')
plt.close()

#P-Pdot plot radio gamma pulsars
plt.figure(4)
plt.scatter(P_radio_gamma,P_dot_radio_gamma,c='purple',marker='o',s=10)
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio gamma pulsars")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_rg.png')
plt.close()

#Positions plot
plt.figure(5)
plt.scatter(x,y,s=2,c='red',alpha=0.4,label='Simulation data',zorder=2)
plt.scatter(y_dirson22,x_dirson22,s=2,c='green',alpha=0.7,label='Dirson et al.(2022)')
#plt.scatter(x2,y2,s=2,c='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.scatter(x_selec,y_selec,s=2,c='blue',alpha=0.5,label='ATNF data') #Canonical pop
plt.scatter([0],[8.5],c='yellow',marker='o',s=20,label='The Sun',zorder=3) #position of the sun
plt.xlim(-30,30)
plt.ylim(-30,30)
#plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend()
plt.savefig('Positions_detected_pulsars.png',dpi=300)
plt.close()

#Distance histogram
plt.figure(6)
plt.hist(distance,bins=20,range=(0,25),edgecolor='black',color='red',alpha=0.5,label='Simulation',zorder=3)
#plt.hist(d2,bins=20,range=(0,25),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.hist(d_selec,bins=20,range=(0,25),edgecolor='black',color='blue',alpha=0.5,label='ATNF data',zorder=2) #Canonical pop
#plt.hist(dist_dirson22,bins=20,range=(0,25),edgecolor='black',color='green',alpha=0.8,label='Dirson et al. (2022)')
plt.legend()
plt.yscale('log')
plt.xlabel('d (kpc)')
plt.ylabel('Frequency')
#plt.title('Distance to earth of the pulsars')
plt.savefig('histo_dist.png',dpi=300)
plt.close()

#Age histogram
plt.figure(7)
plt.hist(log_age,bins=20,range=(1,11),edgecolor='black',color='red',alpha=0.5,label='Simulation')
#plt.hist(log_age2,bins=20,range=(1,11),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.hist(log_age_selec,bins=20,range=(1,11),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Canonical pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
#plt.title('Histogram of the age of the detected pulsars')
plt.savefig('histo_age.png',dpi=300)
plt.close()

#Latitude histogram
plt.figure(8)
plt.hist(latitude,bins=20,range=(-60,60),edgecolor='black',color='red',alpha=0.5,label='Simulation',zorder=3)
plt.hist(lat_dirson22,bins=20,range=(-60,60),edgecolor='black',color='green',alpha=0.5,label='Dirson et al.(2022)',zorder=1)
#plt.hist(latitude2,bins=20,range=(-60,60),edgecolor='black',color='blue',alpha=0.5,label='ATNF data',zorder=2) #Whole pop
plt.hist(latitude_selec,bins=20,range=(-60,60),edgecolor='black',color='blue',alpha=0.5,label='ATNF data',zorder=2) #Canonical pop
plt.yscale('log')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the latitude of the detected pulsars')
plt.savefig('histo_latitude.png',dpi=300)
plt.close()

#Log(P) histogram
plt.figure(9)
plt.hist(log_P,bins=20,range=(-2,1.5),edgecolor='black',color='red',alpha=0.5,label='Simulation')
#plt.hist(log_Pa,bins=20,range=(-2,1.5),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop 
plt.hist(log_P_selec,bins=20,range=(-2,1.5),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Canonical pop
#plt.hist(P_old,bins=20,range=(-2,1.5),edgecolor='black',color='green',alpha=0.5,label='Old pulsars')
plt.legend()
plt.xlabel('Log(P) (P in s)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period of the detected pulsars')
plt.savefig('histo_period.png',dpi=300)
plt.close()

#Log(P_dot) histogram
plt.figure(10)
plt.hist(log_Pdot,bins=20,range=(-20,-10),edgecolor='black',color='red',alpha=0.5,label='Simulation')
#plt.hist(log_Pdot2,bins=20,range=(-20,-10),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.hist(log_Pdot_selec,bins=20,range=(-20,-10),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Canonical pop
plt.legend()
plt.xlabel('Log(Pdot)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period derivative of the detected pulsars')
plt.savefig('histo_pdot.png')
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(11)
plt.hist(cos_alpha,bins=20,range=(-1,1),edgecolor='black',color='red',alpha=0.5,label='cos(alpha)')
plt.hist(cos_alpha0,bins=20,range=(-1,1),edgecolor='black',color='blue',alpha=0.5,label='cos(alpha0)')
plt.legend()
plt.xlabel('cos(alpha) and cos(alpha0)')
plt.ylabel('Frequency')
#plt.title('Histogram of the inclination angle of the detected pulsars')
plt.savefig('histo_cosalpha.png')
plt.close()

#Period of old pulsars : histogram
plt.figure(12)
plt.hist(P_old,bins=20,range=(-2,1.5),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Log(P) (P in s)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period of the detected pulsars with ages betwen 1e7.7 years and 1e8.7')
plt.savefig('histo_period_old.png')
plt.close()

#Relative error between B(P,Pdot) and Bf the decaying magnetic field
plt.figure(13)
plt.hist(err_rel_B,bins=20,range=(0,0.01),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Relative error between the magnetic field computed with decay or with P and Pdot')
plt.ylabel('Frequency')
#plt.title('Histogram of the error on the magnetic field for the detected pulsars')
plt.savefig('histo_err_B.png')
plt.close()

#Velocities
plt.figure(14)
plt.hist(v_norm,bins=20,range=(0,1000),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Velocities of the pulsars in km/s')
plt.ylabel('Frequency')
#plt.title("Histogram of the velocities of the detected pulsars in the simulation")
plt.savefig('histo_velocities.png')
plt.close()

#Comparison between Bfield computed with the decay and with P and Pdot
plt.figure(15)
plt.hist(Bf,bins=20,range=(1e6,1.5e8),edgecolor='black',color='red',alpha=0.5,label='Decaying magnetic field')
plt.hist(B_ppdot,bins=20,range=(1e6,1.5e8),edgecolor='black',color='blue',alpha=0.5,label='Magnetic field computed with P and Pdot')
plt.legend()
plt.xlabel('Magnetic field (B in Tesla)')
plt.ylabel('Frequency')
plt.savefig('histo_Bfield.png')
plt.close()

#Spin-velocity angles histogram
plt.figure(16)
plt.hist(PA,bins=20,range=(0,180),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Frequency')
plt.yscale('log')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle.png',dpi=300)
plt.close()

#Spin-velocity angle (old pulsars) histogram
plt.figure(17)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('spin-velocity angle in degrees for pulsars older than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_old.png')
plt.close()

#Spin-velocity angle (young pulsars) histogram
plt.figure(18)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees for pulsars younger than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_young.png')
plt.close()

#Spin-velocity angle (young pulsars and old) histogram
plt.figure(19)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='black',color='blue',alpha=0.5,label='Age < 10Myr',zorder=2)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='black',color='red',alpha=0.5,label='Age > 10Myr')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_young_and_old.png',dpi=300)
plt.close()

#Plot age=f(spin-vel angle)
plt.figure(20)
plt.scatter(PA,log_age,c='red',marker='o',s=5,label='Simulation data')
plt.xlim(0,180)
plt.ylim(0,12)
#plt.yscale('log')
#plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Log(age) in yr')
plt.legend()
plt.savefig('spin_vel_age_plot.png')
plt.close()

#Initial velocities
plt.figure(21)
plt.hist(v0_norm,bins=20,range=(0,1000),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Birth kick velocities of the pulsars in km/s')
plt.ylabel('Frequency')
#plt.title("Histogram of the velocities of the detected pulsars in the simulation")
plt.savefig('histo_BK_velocities.png')
plt.close()

#with open("save_coord.txt","r") as f:
#    data_x_y=re.findall(reg_4,f.read())

#x0_pulsar,y0_pulsar=[],[]
#for i in range(int(len(data_x_y)/2)):
#    x0_pulsar.append(float(data_x_y[2*i]))
#    y0_pulsar.append(float(data_x_y[2*i+1]))

#plt.figure(22)
#plt.xlim(-20,20)
#plt.ylim(-20,20)
#plt.scatter(x0_pulsar,y0_pulsar,c='blue',marker='*',s=0.0000001,label='positions of pulsars')
#plt.savefig("initial_positions_pulsars.png",dpi=300)
#plt.close()

#Death line alone
#plt.figure(22)
#plt.plot(P_line,P_dot_line_CR93)
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig("death_line.png")
#plt.close()

#Binit histogram
#plt.figure(12)
#plt.hist(B_init,bins=150,range=(1e7,9e8),edgecolor='black',color='red',alpha=0.5,label='Simulation')
#plt.legend()
#plt.xlabel('Inital magnetic field B')
#plt.ylabel('Frequency')
#plt.title('Histogram of the initial magnetic field of the pulsars')
#plt.savefig('B_init.png')
