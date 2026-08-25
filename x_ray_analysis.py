import numpy as np
import matplotlib.pyplot as plt
import re
from mpl_toolkits import mplot3d
from scipy.stats import kstest,ks_2samp,linregress
from scipy.stats import mannwhitneyu
from matplotlib.colors import LogNorm
from astropy.table import Table
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC

#Get data from X-ray catalog (no BB column)
df=pd.read_excel('X_ray_data_wu_et_al.ods')
#df["Pdot"]=(df["Pdot"].str.replace("\u00D710","E",regex=False).astype(float))
#df["LX"]=(df["LX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LSX"]=(df["LSX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LHX"]=(df["LHX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LG"]=(df["LG"].str.replace("\u00D710","E",regex=False).astype(str))
data_X=Table.from_pandas(df)
for i in range(len(data_X['Sorting_distance'])):
    if str(data_X['Sorting_distance'][i])[0]!='<':
        data_X['Sorting_distance'][i]=float(data_X['Sorting_distance'][i])
    else:
        val=str(data_X['Sorting_distance'][i])[1:]
        data_X['Sorting_distance'][i]=val
        data_X['Sorting_distance'][i]=float(data_X['Sorting_distance'][i])
mask=data_X['Sorting_distance'] < 25
data_X_filtered=data_X[mask]
data_X=data_X_filtered
data_X['LX_upper_limit'] = [str(x).startswith('<') for x in data_X['LX']]
data_X['LX'] = [float(str(x).replace('<','')) for x in data_X['LX']]
data_X['LX']=data_X['LX']*1e-7
#df.to_excel("X_ray_data_wu_et_al.ods", engine="odf", index=False)
#print(data_X)

#Get data from X-ray catalog (BB column)
df2=pd.read_excel('X_ray_data_wu_et_al_colspec.ods')
#df["Pdot"]=(df["Pdot"].str.replace("\u00D710","E",regex=False).astype(float))
#df["LX"]=(df["LX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LSX"]=(df["LSX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LHX"]=(df["LHX"].str.replace("\u00D710","E",regex=False).astype(str))
#df["LG"]=(df["LG"].str.replace("\u00D710","E",regex=False).astype(str))
data_X2=Table.from_pandas(df2)
for i in range(len(data_X2['Sorting_distance'])):
    if str(data_X2['Sorting_distance'][i])[0]!='<':
        data_X2['Sorting_distance'][i]=float(data_X2['Sorting_distance'][i])
    else:
        val=str(data_X2['Sorting_distance'][i])[1:]
        data_X2['Sorting_distance'][i]=val
        data_X2['Sorting_distance'][i]=float(data_X2['Sorting_distance'][i])
mask=data_X2['Sorting_distance'] < 25
data_X2_filtered=data_X2[mask]
data_X2=data_X2_filtered
data_X2['LX_upper_limit'] = [str(x).startswith('<') for x in data_X2['LX']]
data_X2['LX'] = [float(str(x).replace('<','')) for x in data_X2['LX']]
data_X2['LX']=data_X2['LX']*1e-7
mask=data_X2['Spectrum_BB'] == 1
data_X2_fil=data_X2[mask]
data_X2=data_X2_fil
#df.to_excel("X_ray_data_wu_et_al.ods", engine="odf", index=False)
#print(data_X2)

#Count the pulsars (observed) -> without filter on BB spectrum component 
#nb_X,nb_RX,nb_GX,nb_RGX,nb_pulse=0,0,0,0,0
#for i in range(len(data_X["P"])):
#    if data_X['X_pulsation'][i]==1:
#        nb_pulse+=1
#    if data_X["Type"][i]==4:
#        nb_X+=1
#    elif data_X["Type"][i]==1:
#        nb_RX+=1
#    elif data_X["Type"][i]==2:
#        nb_GX+=1
#    elif data_X["Type"][i]==3:
#        nb_RGX+=1
#print(f'----OBSERVATIONS----')
#print(f'Number of X-ray only pulsars: {nb_X}\nNumber of Radio/X-ray pulsars: {nb_RX}\nNumber of gamma-ray/X-ray pulsars: {nb_GX}\nNumber of Radio/Gamma-ray/X-ray pulsars: {nb_RGX}')
#print(f'Number of pulsating X-ray sources: {nb_pulse}')

#with open("info_supp_obs.txt","a") as f:
#    f.write(f'{nb_GX+nb_RGX}\n{nb_pulse}\n')

#Count the pulsars (observed) -> with filter on BB spectrum component
nb_X2,nb_RX2,nb_GX2,nb_RGX2,nb_pulse2=0,0,0,0,0
for i in range(len(data_X2["P"])):
    if data_X2['X_pulsation'][i]==1:
        nb_pulse2+=1
    if data_X2["Type"][i]==4:
        nb_X2+=1
    elif data_X2["Type"][i]==1:
        nb_RX2+=1
    elif data_X2["Type"][i]==2:
        nb_GX2+=1
    elif data_X2["Type"][i]==3:
        nb_RGX2+=1
print(f'----OBSERVATIONS----')
print(f'Number of X-ray only pulsars: {nb_X2}\nNumber of Radio/X-ray pulsars: {nb_RX2}\nNumber of gamma-ray/X-ray pulsars: {nb_GX2}\nNumber of Radio/Gamma-ray/X-ray pulsars: {nb_RGX2}')
print(f'Number of pulsating X-ray sources: {nb_pulse2}')

with open("info_supp_obs.txt","a") as f:
    f.write(f'{nb_GX2+nb_RGX2}\n{nb_pulse2}\n')

#Simulation data
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")

with open("x_file.txt","r") as f:
    data=re.findall(reg_1,f.read())

with open("x_file.txt","r") as f:
    data_type=re.findall(reg_2,f.read())

with open("x_file2.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("PF_check.txt","r") as f:
    data_PF=re.findall(reg_3,f.read())

Fxmax,Fxmin,PF=[],[],[]
for i in range(int(len(data_PF)/3)):
    PF.append(float(data_PF[3*i+0]))
    Fxmax.append(float(data_PF[3*i+1]))
    Fxmin.append(float(data_PF[3*i+2]))

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

count_x_sim,count_rgx_sim,count_rx_sim,count_gx_sim=0,0,0,0
for i in range(len(type_pulsar)):
    if (type_pulsar[i]==1):
        count_rx_sim+=1
    elif (type_pulsar[i]==2):
        count_gx_sim+=1
    elif (type_pulsar[i]==3):
        count_rgx_sim+=1
    elif (type_pulsar[i]==4):
        count_x_sim+=1

print(f'----SIMULATIONS----')
print(f'Number of X-ray only pulsars: {count_x_sim}\nNumber of Radio/X-ray pulsars: {count_rx_sim}\nNumber of gamma-ray/X-ray pulsars: {count_gx_sim}\nNumber of Radio/Gamma-ray/X-ray pulsars: {count_rgx_sim}')

age_year=[]
for i in range(len(age)):
    age_year.append(age[i]/(365.25*24*60*60))

alpha=[]
for i in range(len(cos_alpha)):
    alpha.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))

T_n,T_s,r_hn,cos_in,F_x,xi,r_hs=[],[],[],[],[],[],[]
cos_is,theta_n,phi_n,theta_s,phi_s,mu_hot_ang1,mu_hot_ang2=[],[],[],[],[],[],[]
for i in range(int(len(data2)/14)):
    cos_in.append(float(data2[14*i])) #cos of the angle between magnetic axis and the line of sight (north)
    T_n.append(float(data2[14*i+1])) #Temperature of the hot spot in K
    r_hn.append(float(data2[14*i+2])) #Radius of the hot spot in meter
    F_x.append(float(data2[14*i+3])) #Flux in thermal X-ray in W.m^-2
    xi.append(float(data2[14*i+4])*180/np.pi) #Viewing angle in degree
    cos_is.append(float(data2[14*i+5])) #cos of the angle between magnetic axis and the line of sight (south)
    theta_n.append(float(data2[14*i+6])) #Position in spherical coordinate in rad (theta and north)
    phi_n.append(float(data2[14*i+7])) #Position in spherical coordinate in rad (phi and north)
    theta_s.append(float(data2[14*i+8])) #Position in spherical coordinate in rad (theta and south)
    phi_s.append(float(data2[14*i+9])) #Position in spherical coordinate in rad (phi and south)
    mu_hot_ang1.append(float(data2[14*i+10])*180/np.pi) #Angle in degree between magnetic axis and hotspot (north)
    mu_hot_ang2.append(float(data2[14*i+11])*180/np.pi) #Angle in degree between magnetic axis and hotspot (south)
    r_hs.append(float(data2[14*i+12])) #Radius of the hot spot in meter (south)
    T_s.append(float(data2[14*i+13])) #Temperature of the hot spot in K (south)

cos_in=[None if x==1e55 else x for x in cos_in]
cos_is=[None if x==1e55 else x for x in cos_is]
theta_n=[None if x==1e55 else x for x in theta_n]
theta_s=[None if x==1e55 else x for x in theta_s]
phi_n=[None if x==1e55 else x for x in phi_n]
phi_s=[None if x==1e55 else x for x in phi_s]
mu_hot_ang1=[None if x==1e55 else x for x in mu_hot_ang1]
mu_hot_ang2=[None if x==1e55 else x for x in mu_hot_ang2]
r_hn=[None if x==1e55 else x for x in r_hn]
r_hs=[None if x==1e55 else x for x in r_hs]
T_n=[None if x==1e55 else x for x in T_n]
T_s=[None if x==1e55 else x for x in T_s]

all_rh=r_hn+r_hs
all_theta=theta_n+theta_s
all_phi=phi_n+phi_s
all_mu_hot_ang=mu_hot_ang1+mu_hot_ang2
all_cos_i=cos_in+cos_is
all_T=T_n+T_s

all_theta=np.array(all_theta,dtype=float)
all_phi=np.array(all_phi,dtype=float)
all_cos_i = np.array(all_cos_i, dtype=float)
all_rh=np.array(all_rh,dtype=float)
all_T=np.array(all_T,dtype=float)

count_big_T=(all_T>7e6).sum()
print(count_big_T)

Lx_BB=[]
sigma=5.67e-8
for i in range(len(F_x)):
    if (r_hn[i]!=None and r_hs[i]!=None):
        Lx_BB.append((np.pi*(r_hn[i])**2*sigma*T_n[i]**4)+(np.pi*(r_hs[i])**2*sigma*T_s[i]**4))
    elif (r_hs[i]!=None and r_hn[i]==None):
        Lx_BB.append(np.pi*(r_hs[i])**2*sigma*T_s[i]**4)
    elif (r_hs[i]==None and r_hn[i]!=None):
        Lx_BB.append(np.pi*(r_hn[i])**2*sigma*T_n[i]**4)

Lx_abs_redshift=[]
for i in range(len(F_x)):
    Lx_abs_redshift.append(F_x[i]*4*np.pi*(distance[i]*3.086e19)**2)

#Plot the death line of the article from Mitra et al. (2019)
Inertia=1e38
R_NS=12000
mu_0=1.25663706212e-6
c_light=2.997924858e8
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
P_line=[10**(np.log10(i)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line=[10**np.log10(const*(i**2)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line4=[Pdot_line[i]*10**(-0.55) for i in range(len(Pdot_line))]
Pdot_line5=[Pdot_line[i]*10**(1.15) for i in range(len(Pdot_line))]
Pdot_line2=[10**np.log10(const_min*(i**2)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line3=[10**np.log10(const_max*(i**2)) for i in np.arange(1e-2,3e1,0.01)]

#Plot the line with the critical magnetic field 4.4e9 T, distinguishing magnetar from canonical pulsars
B_crit=4.4e9
B0=3.2e15 #constant to compute the surface magnetic field of the observations
B_linecrit=[10**np.log10(((B_crit/B0)**2)*(i**(-1))) for i in np.arange(1e-2,3e1,0.01)]

#Plot the Edot lines
P_dot_Edot1e22W=[(1e22*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e23W=[(1e23*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e25W=[(1e25*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e28W=[(1e28*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e31W=[(1e31*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]

#Plot the B lines
P_dot_B1e4=[(1e4/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e5=[(1e5/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e6=[(1e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B5e6=[(5e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e7=[(1e7/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e8=[(1e8/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e9=[(1e9/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B4_4e9=[(4.4e9/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]

#Edot
Edot_obs,Edot_sim=[],[]
for i in range(len(P)):
    Edot_sim.append(4*np.pi**2*Inertia*P_dot[i]*P[i]**(-3))

for i in range(len(data_X2['P'])):
    Edot_obs.append(4*np.pi**2*Inertia*data_X2['Pdot'][i]*data_X2['P'][i]**(-3))

#KS test 1D python (all the data)
#KS_test=kstest(P_dot,data_X['Pdot'])
#test_stat_all1=KS_test.statistic
#p_value_all1=KS_test.pvalue
#print("----ALL THE DATA----")
#print(f"d_value of Pdot KS test = {test_stat_all1}")
#print(f"p_value of Pdot={p_value_all1}")

#KS_test=kstest(P,data_X['P'])
#test_stat_all2=KS_test.statistic
#p_value_all2=KS_test.pvalue
#print(f"d_value of P KS test = {test_stat_all2}")
#print(f"p_value of P={p_value_all2}")

condition = [Pdot2 < Pdot3 for Pdot2, Pdot3 in zip(Pdot_line2,Pdot_line3)]

#P-Pdot plot all pulsars
plt.figure(1)
plt.scatter(P,P_dot,c='red',marker='o',s=5,label='Simulation data',zorder=2)
plt.scatter(data_X2['P'],data_X2['Pdot'],c='blue',marker='o',s=5,label='X-ray catalog of Xu et al. (2025)',zorder=1) #Only canonical pop 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=1) #Death line Mitra et al. 2019
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley' )

#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e22W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1800]*0.85,P_dot_Edot1e22W[1800]*1.05,r'$\dot{E} = 10^{22}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[800]*0.85,P_dot_Edot1e25W[800]*1.05,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[77]*0.85,P_dot_Edot1e28W[77]*1.05,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[7]*0.85,P_dot_Edot1e31W[7]*1.05,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[80]*0.6,P_dot_B1e6[80]*0.95,r'$B =10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[1500]*0.6,P_dot_B5e6[1500]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[2000]*0.6,P_dot_B1e7[2000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[2300]*0.6,P_dot_B1e8[2300]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.plot(P_line,P_dot_B1e9,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[2500]*0.6,P_dot_B1e9[2500]*0.95,r'$B= 10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B4_4e9,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[2500]*0.6,P_dot_B4_4e9[2500]*0.95,r'$B= 4.4\times10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)

plt.xlim(1e-2,3e1)
plt.ylim(1e-19,1e-11)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(loc='lower left',fontsize='x-small')
plt.savefig('P_Pdot_plot_X.pdf',dpi=300)
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(2)
plt.hist(cos_alpha,bins=20,range=(-1,1),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'cos($\alpha$)')
plt.hist(cos_alpha0,bins=20,range=(-1,1),edgecolor='blue',color='blue',histtype='step',alpha=0.5,label=r'cos($\alpha_0$)')
plt.legend()
plt.xlabel(r'cos($\alpha$) and cos($\alpha_0$)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_cosalphaX.pdf',dpi=300)
plt.close()

#cos(i) histogram
#plt.figure(3)
#plt.hist(cos_i,bins=20,range=(-1,1),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'cos($i$)')
#plt.legend()
#plt.xlabel(r'cos($i$)')
#plt.ylabel('Number of pulsars')
#plt.savefig('histo_cos_i.pdf',dpi=300)
#plt.close()

#Lx histogram
plt.figure(4)
plt.hist(np.log10(data_X2['LX']),bins=20,range=(21,28),edgecolor='blue',color='blue',histtype='step',alpha=0.5,label=r'X-ray catalog of Xu et al. (2025)')
plt.hist(np.log10(Lx_abs_redshift),bins=20,range=(21,28),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'$L_X$ in logscale (W)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_LX.pdf',dpi=300)
plt.close()

#T histogram
plt.figure(5)
plt.hist(np.log10(all_T),bins=20,edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
#plt.xscale('log')
plt.xlabel(r'Temperature of the hot spot in logscale (K)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_T.pdf',dpi=300)
plt.close()

#r_h histogram
plt.figure(6)
plt.hist(all_rh,bins=20,edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'Radius of the hot spot (m)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_rh.pdf',dpi=300)
plt.close()

#Zeta=f(chi) all
plt.figure(7)
plt.scatter(alpha,xi,c='red',marker='o',s=5,label='Simulation')
plt.ylabel(r'$\zeta$ in °')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.legend(fontsize='small')
plt.savefig('zeta_chi_X.pdf',dpi=300)
plt.close()

#F_x histogram
plt.figure(8)
plt.hist(np.log10(F_x),bins=20,range=(-20,-10),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'X-ray flux in logscale (W.m$^{-2}$)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_flux.pdf',dpi=300)
plt.close()

#Linear regression Lx=f(Edot) (obs)
a_lx,b_lx,r_value,p_value,std_err=linregress(np.log10(Edot_obs),np.log10(data_X2['LX']))
reglinx = np.logspace(np.log10(min(Edot_obs)), np.log10(max(Edot_obs)), num=len(Edot_obs))
regliny = 10**(a_lx * np.log10(reglinx) + b_lx)
residuals_y= np.log10(data_X2['LX']) - np.log10(regliny)
sigma_y = np.sqrt(np.sum(residuals_y**2) / (len(Edot_obs) - 2))
Sxx=np.sum((np.log10(Edot_obs) - np.mean(np.log10(Edot_obs)))**2)
std_err_b=sigma_y*np.sqrt(1.0/len(Edot_obs)+((np.mean(np.log10(Edot_obs)))**2/Sxx))
a_lxplus=a_lx+std_err
b_lxplus=b_lx+std_err_b
a_lxminus=a_lx-std_err
b_lxminus=b_lx-std_err_b
reglinyplus = 10**(a_lx * np.log10(reglinx) + b_lxplus)
reglinyminus = 10**(a_lx * np.log10(reglinx) + b_lxminus)

print(a_lxplus)
print(a_lxminus)

#Linear regression Lx=f(Edot) (sim: abs+redshift)
if (len(Edot_sim)>=2 and len(Lx_abs_redshift)>=2):
    a_lx2,b_lx2,r_value2,p_value2,std_err2=linregress(np.log10(Edot_sim),np.log10(Lx_abs_redshift))
    reglinx2 = np.logspace(np.log10(min(Edot_sim)), np.log10(max(Edot_sim)), num=len(Edot_sim))
    regliny2 = 10**(a_lx2 * np.log10(reglinx2) + b_lx2)
    residuals_y2= np.log10(Lx_abs_redshift) - np.log10(regliny2)
    sigma_y2 = np.sqrt(np.sum(residuals_y2**2) / (len(Edot_sim) - 2))
    Sxx2=np.sum((np.log10(Edot_sim) - np.mean(np.log10(Edot_sim)))**2)
    std_err_b2=sigma_y2*np.sqrt(1.0/len(Edot_sim)+((np.mean(np.log10(Edot_sim)))**2/Sxx2))
    a_lxplus2=a_lx2+std_err2
    b_lxplus2=b_lx2+std_err_b2
    a_lxminus2=a_lx2-std_err2
    b_lxminus2=b_lx2-std_err_b2
    reglinyplus2 = 10**(a_lx2 * np.log10(reglinx2) + b_lxplus2)
    reglinyminus2 = 10**(a_lx2 * np.log10(reglinx2) + b_lxminus2)
else:
    a_lx2,b_lx2=0,0

#Linear regression Lx=f(Edot) (sim: bolometric)
if (len(Edot_sim)>=2 and len(Lx_BB)>=2):
    a_lx3,b_lx3,r_value3,p_value3,std_err3=linregress(np.log10(Edot_sim),np.log10(Lx_BB))
    reglinx3 = np.logspace(np.log10(min(Edot_sim)), np.log10(max(Edot_sim)), num=len(Edot_sim))
    regliny3 = 10**(a_lx3 * np.log10(reglinx3) + b_lx3)
    residuals_y3= np.log10(Lx_BB) - np.log10(regliny3)
    sigma_y3 = np.sqrt(np.sum(residuals_y3**2) / (len(Edot_sim) - 2))
    Sxx3=np.sum((np.log10(Edot_sim) - np.mean(np.log10(Edot_sim)))**2)
    std_err_b3=sigma_y3*np.sqrt(1.0/len(Edot_sim)+((np.mean(np.log10(Edot_sim)))**2/Sxx3))
    a_lxplus3=a_lx3+std_err3
    b_lxplus3=b_lx3+std_err_b3
    a_lxminus3=a_lx3-std_err3
    b_lxminus3=b_lx3-std_err_b3
    reglinyplus3 = 10**(a_lx3 * np.log10(reglinx3) + b_lxplus3)
    reglinyminus3 = 10**(a_lx3 * np.log10(reglinx3) + b_lxminus3)
else:
    a_lx3,b_lx3=0,0

with open("info_supp.txt", "a") as f:
    if (len(Edot_sim)>=2 and len(Lx_BB)>=2):
        f.write(f'{a_lx3}\n')
    else: 
        f.write('0\n')

with open("info_supp_obs.txt","a") as f:
    f.write(f'{a_lx}\n')

#Lx=f(Edot)
plt.figure(9)
plt.scatter(Edot_sim,Lx_BB,c='green',s=5,marker='o',label=r'Sim, bolometric $L_X$')
plt.scatter(Edot_sim,Lx_abs_redshift,c='red',s=5,marker='o',label='Simulation')
plt.scatter(Edot_obs,data_X2['LX'],c='blue',s=5,marker='o',label='X-ray catalog of Xu et al. (2025)')
plt.plot(reglinx,regliny,linestyle='-',label=r'$\log$($L_X$) = %.2f $\log(\dot{E})$ + %.2f (obs)' % (a_lx, b_lx),c='blue')
if (len(Edot_sim)>=2 and len(Lx_BB)>=2):
    plt.plot(reglinx2,regliny2,linestyle='-',label=r'$\log$($L_X$) = %.2f $\log(\dot{E})$ + %.2f (sim)' % (a_lx2, b_lx2),c='red')
    plt.plot(reglinx3,regliny3,linestyle='-',label=r'$\log$($L_X$) = %.2f $\log(\dot{E})$ + %.2f (sim bolo)' % (a_lx3, b_lx3),c='green')
    plt.fill_between(reglinx2,reglinyminus2, reglinyplus2, color='red', alpha=0.1)
plt.fill_between(reglinx, reglinyminus, reglinyplus, color='blue', alpha=0.1)
plt.legend(fontsize='small')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$L_X$ (W)')
plt.xlabel(r'$\dot{E}$ in W')
plt.xlim(1e23,1e33)
plt.ylim(1e20,1e30)
plt.savefig('Lx_Edot.pdf',dpi=300)
plt.close()

#age histogram
plt.figure(10)
plt.hist(np.log10(age_year),bins=20,range=(0,9),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'Age in log space (yr)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_ageX.pdf',dpi=300)
plt.close()

#Position of the hotspot
plt.figure(11)
ax=plt.axes(projection='3d')
#Sphere
us = np.linspace(0, 2*np.pi, 100)
vs = np.linspace(0, np.pi, 100)
xsph = np.outer(np.cos(us), np.sin(vs))
ysph = np.outer(np.sin(us), np.sin(vs))
zsph = np.outer(np.ones_like(us), np.cos(vs))
ax.plot_surface(xsph, ysph, zsph, alpha=0.1)
#Cartesian coordinates
xhot = np.sin(all_theta) * np.cos(all_phi)
yhot = np.sin(all_theta) * np.sin(all_phi)
zhot = np.cos(all_theta)
ax.scatter(xhot, yhot, zhot, s=5)
ax.set_box_aspect([1,1,1])
plt.savefig('loc_hotspot_on_sphere.pdf',dpi=300)
plt.close()

#Angle between magnetic axis and hotspot histogram
plt.figure(12)
plt.hist(all_mu_hot_ang,bins=20,range=(0,180),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'Angle between the magnetic axis and the hotspots centers (°)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_mu_hot_ang.pdf',dpi=300)
plt.close()

#Angle between magnetic axis and hotspot histogram
plt.figure(13)
plt.hist(all_cos_i,bins=20,range=(-1,1),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'$\cos(i)$')
plt.ylabel('Number of pulsars')
plt.savefig('histo_all_cosi.pdf',dpi=300)
plt.close()

#Position of the hotspot 2D projection
all_theta_proj=np.pi/2-all_theta
all_phi_proj=all_phi
all_phi_proj[all_phi_proj>np.pi]-=2*np.pi
plt.figure(14)
plt.axes(projection="mollweide")
plt.scatter(all_phi_proj,all_theta_proj,s=10,c='red',alpha=1,label='Simulation data',zorder=2)
plt.grid(alpha=0.3)
plt.legend(loc="upper right", bbox_to_anchor=(1.05, 1.15), borderaxespad=0.)
plt.xlabel('Longitude in °')
plt.ylabel('Colatitude in °')
plt.savefig('loc_hotspot_2Dproj.pdf',dpi=300)
plt.close()

#PF histogram
plt.figure(15)
plt.hist(PF,bins=20,range=(0,1),edgecolor='red',color='red',histtype='step',alpha=0.5,label=r'Simulation')
plt.legend()
plt.xlabel(r'Pulsed fraction')
plt.ylabel('Number of pulsars')
plt.savefig('histo_PF.pdf',dpi=300)
plt.close()
