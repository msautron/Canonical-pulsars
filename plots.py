import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d
from scipy.stats import kstest,ks_2samp
from scipy.stats import mannwhitneyu
from matplotlib.colors import LogNorm
from astropy.table import Table
import pandas as pd

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out 
Pa_gamma,P_dota_gamma,da_gamma,za_gamma,xa_gamma,ya_gamma,agea_gamma,E_dota_gamma=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue of the gamma pulsars only
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
test_l=[]
Ba,B2,Ba_gamma=[],[],[]
B0=3.2e15 #constant to compute the surface magnetic field of the observations
P_selected,Pdot_selected=[],[]
#B_init=[]

#Put the data at the right place
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("[-+]?\d*[.]\d*[Ee]*[-+]*\d*")
reg_4=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")

with open("ATNF_canonical_pulsars_gamma.txt","r") as f:
    data_only_gamma=re.findall(reg_3,f.read())

with open("ATNF_canonical_pulsars_radio.txt","r") as f:
    data_only_r=re.findall(reg_3,f.read())

with open("ATNF_canonical_pulsars.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("ATNF_canonical_pulsars_rg.txt","r") as f:
    data_only_rg=re.findall(reg_3,f.read())

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

with open("Fg_flux.txt","r") as f:
    fg_data=re.findall(reg_4,f.read())

with open("gamma_peak_sep.txt","r") as f:
    g_peak_sep_data=re.findall(reg_4,f.read())

with open("nb_orbit.txt","r") as f:
    nb_orbit_data=re.findall(reg_4,f.read())

with open("lat_dist_noISM.txt","r") as f :
    lines = f.readlines()

df=pd.read_excel('/home/matteo.sautron/Documents/cuda/canonical_pulsars/pulsar_population_repo/3PC_Catalog_20230803.xls')
data_3PC=Table.from_pandas(df)

#Get the flux from the 3PC catalog for the canonical pulsars
flux_3PC_cano,g_peak_sep_obs=[],[]
for i in range(len(data_3PC['G100'])):
    B_3PC=3.2e15*np.sqrt(data_3PC['P0'][i]*data_3PC['P1'][i])
    if B_3PC<=4.4e9 and B_3PC>=6e6 and data_3PC['G100'][i]!='*':
        flux_3PC_cano.append(float(data_3PC['G100'][i])*1e-3)
    if B_3PC<=4.4e9 and B_3PC>=6e6 and data_3PC['PKSEP'][i]!='*':
        if data_3PC['PKSEP'][i] <= 0.5: 
            g_peak_sep_obs.append(float(data_3PC['PKSEP'][i]))
        else:
            g_peak_sep_obs.append(1.0-float(data_3PC['PKSEP'][i]))

print(f'Number of canonical gamma pulsars : {len(g_peak_sep_obs)}\n')

#Get the data of distance and latitude when the ISM is not taken into account
lat_noISM,dist_noISM,x_noISM,y_noISM=[],[],[],[]
for i,line in enumerate(lines):
    values=line.strip().split()
    lat_noISM.append(float(values[0]))
    dist_noISM.append(float(values[1]))
    x_noISM.append(float(values[2]))
    y_noISM.append(float(values[3]))

#Get the data from Dirson et al. (2022) for xy coord and latitude
x_dirson22,y_dirson22,lat_dirson22,dist_dirson22=[],[],[],[]
for i in range(int(len(data1_dirson22)/2)):
    x_dirson22.append(float(data1_dirson22[2*i]))
    y_dirson22.append(float(data1_dirson22[2*i+1]))

for i in range(int(len(data2_dirson22))):
    lat_dirson22.append(float(data2_dirson22[i]))

for i in range(int(len(data3_dirson22))):
    dist_dirson22.append(float(data3_dirson22[i]))

#Get the data on the Fg flux, allowing then to plot the new pulsars detected in gamma with a sensitivity 10 times greater than before
Fg_flux=[]
for i in range(int(len(fg_data))):
    Fg_flux.append(float(fg_data[i]))

#Get the data about the gamma-ray peak separation
g_peak_sep_sim=[]
for i in range(int(len(g_peak_sep_data))):
    g_peak_sep_sim.append(float(g_peak_sep_data[i]))

#Get the info about the nb of orbits made by each NS (approximately) 
nb_orbit=[]
for i in range(int(len(nb_orbit_data))):
    nb_orbit.append(float(nb_orbit_data[i]))

#Get the data of the ATNF catalogue
for i in range(int(len(data2)/9)):
    Pa+=[float(data2[i*9])] #Period of the rotation of the pulsar in seconds
    P_dota+=[float(data2[i*9+1])] #Period derivative of the rotation of the pulsar no units
    da+=[float(data2[i*9+2])] #Distance to us in kpc
    za+=[float(data2[i*9+3])] #Z position in the galactocentric frame in kpc
    xa+=[float(data2[i*9+4])] #X position in the galactocentric frame in kpc
    ya+=[float(data2[i*9+5])] #age of the pulsar in years
    agea+=[float(data2[i*9+6])] #age of the pulsars in seconds
    E_dota+=[float(data2[i*9+7])]  #Spin down power of the pulsar in ergs/s
    Ba+=[float(data2[i*9+8])] #surface magnetic field in T

#Get the data of the ATNF catalogue of the gamma pulsars only
for i in range(int(len(data_only_gamma)/9)):
    Pa_gamma+=[float(data_only_gamma[i*9])] #Period of the rotation of the pulsar in seconds
    P_dota_gamma+=[float(data_only_gamma[i*9+1])] #Period derivative of the rotation of the pulsar no units
    da_gamma+=[float(data_only_gamma[i*9+2])] #Distance to us in kpc
    za_gamma+=[float(data_only_gamma[i*9+3])] #Z position in the galactocentric frame in kpc
    xa_gamma+=[float(data_only_gamma[i*9+4])] #X position in the galactocentric frame in kpc
    ya_gamma+=[float(data_only_gamma[i*9+5])] #age of the pulsar in years
    agea_gamma+=[float(data_only_gamma[i*9+6])] #age of the pulsars in seconds
    E_dota_gamma+=[float(data_only_gamma[i*9+7])]  #Spin down power of the pulsar in ergs/s
    Ba_gamma+=[float(data_only_gamma[i*9+8])] #surface magnetic field in T

#Get the data of the ATNF catalogue of the radio-gamma pulsars only
Pa_rg,P_dota_rg,da_rg,za_rg,xa_rg,ya_rg,agea_rg,E_dota_rg,Ba_rg=[],[],[],[],[],[],[],[],[]
for i in range(int(len(data_only_rg)/9)):
    Pa_rg+=[float(data_only_rg[i*9])] #Period of the rotation of the pulsar in seconds
    P_dota_rg+=[float(data_only_rg[i*9+1])] #Period derivative of the rotation of the pulsar no units
    da_rg+=[float(data_only_rg[i*9+2])] #Distance to us in kpc
    za_rg+=[float(data_only_rg[i*9+3])] #Z position in the galactocentric frame in kpc
    xa_rg+=[float(data_only_rg[i*9+4])] #X position in the galactocentric frame in kpc
    ya_rg+=[float(data_only_rg[i*9+5])] #age of the pulsar in years
    agea_rg+=[float(data_only_rg[i*9+6])] #age of the pulsars in seconds
    E_dota_rg+=[float(data_only_rg[i*9+7])]  #Spin down power of the pulsar in ergs/s
    Ba_rg+=[float(data_only_rg[i*9+8])] #surface magnetic field in T

#Get the data of the ATNF catalogue of the radio pulsars only
Pa_r,P_dota_r,da_r,za_r,xa_r,ya_r,agea_r,E_dota_r,Ba_r=[],[],[],[],[],[],[],[],[]
for i in range(int(len(data_only_r)/9)):
    Pa_r+=[float(data_only_r[i*9])] #Period of the rotation of the pulsar in seconds
    P_dota_r+=[float(data_only_r[i*9+1])] #Period derivative of the rotation of the pulsar no units
    da_r+=[float(data_only_r[i*9+2])] #Distance to us in kpc
    za_r+=[float(data_only_r[i*9+3])] #Z position in the galactocentric frame in kpc
    xa_r+=[float(data_only_r[i*9+4])] #X position in the galactocentric frame in kpc
    ya_r+=[float(data_only_r[i*9+5])] #age of the pulsar in years
    agea_r+=[float(data_only_r[i*9+6])] #age of the pulsars in seconds
    E_dota_r+=[float(data_only_r[i*9+7])]  #Spin down power of the pulsar in ergs/s
    Ba_r+=[float(data_only_r[i*9+8])] #surface magnetic field in T

#Create lists with the data of all gamma pulsars (gamma + radio-gamma)
Pa_gall,P_dota_gall,da_gall,za_gall,xa_gall,ya_gall,agea_gall,E_dota_gall,Ba_gall=[],[],[],[],[],[],[],[],[]
Pa_gall=Pa_rg+Pa_gamma
P_dota_gall=P_dota_rg+P_dota_gamma
da_gall=da_rg+da_gamma
za_gall=za_rg+za_gamma
xa_gall=xa_rg+xa_gamma
ya_gall=ya_rg+ya_gamma
agea_gall=agea_rg+agea_gamma
E_dota_gall=E_dota_rg+E_dota_gamma
Ba_gall=Ba_rg+Ba_gamma

Pa_gall_no5,P_dota_gall_no5=[],[]
for i in range(len(Pa_gall)):
    if (Pa_gall[i]<5):
        Pa_gall_no5+=[Pa_gall[i]]
        P_dota_gall_no5+=[P_dota_gall[i]]

#Data we are really using, getting rid of lines where we miss a single info
P2_g,P_dot2_g,d2_g,z2_g,x2_g,y2_g,age2_g,E_dot2_g,B2_g=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dota_gamma)):
    if P_dota_gamma[i]!=0.0 and da_gamma[i]!=0.0 and E_dota_gamma[i]!=0 and da_gamma[i]<=25:
        P2_g+=[Pa_gamma[i]]
        P_dot2_g+=[P_dota_gamma[i]]
        d2_g+=[da_gamma[i]]
        z2_g+=[za_gamma[i]]
        x2_g+=[xa_gamma[i]]
        y2_g+=[ya_gamma[i]]
        age2_g+=[agea_gamma[i]]
        E_dot2_g+=[E_dota_gamma[i]]
        B2_g+=[Ba_gamma[i]]

P2_rg,P_dot2_rg,d2_rg,z2_rg,x2_rg,y2_rg,age2_rg,E_dot2_rg,B2_rg=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dota_rg)):
    if P_dota_rg[i]!=0.0 and da_rg[i]!=0.0 and E_dota_rg[i]!=0 and da_rg[i]<=25:
        P2_rg+=[Pa_rg[i]]
        P_dot2_rg+=[P_dota_rg[i]]
        d2_rg+=[da_rg[i]]
        z2_rg+=[za_rg[i]]
        x2_rg+=[xa_rg[i]]
        y2_rg+=[ya_rg[i]]
        age2_rg+=[agea_rg[i]]
        E_dot2_rg+=[E_dota_rg[i]]
        B2_rg+=[Ba_rg[i]]

P2_gall,P_dot2_gall,d2_gall,z2_gall,x2_gall,y2_gall,age2_gall,E_dot2_gall,B2_gall=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dota_gall)):
    if P_dota_gall[i]!=0.0 and da_gall[i]!=0.0 and E_dota_gall[i]!=0 and da_gall[i]<=25:
        P2_gall+=[Pa_gall[i]]
        P_dot2_gall+=[P_dota_gall[i]]
        d2_gall+=[da_gall[i]]
        z2_gall+=[za_gall[i]]
        x2_gall+=[xa_gall[i]]
        y2_gall+=[ya_gall[i]]
        age2_gall+=[agea_gall[i]]
        E_dot2_gall+=[E_dota_gall[i]]
        B2_gall+=[Ba_gall[i]]

P2_r,P_dot2_r,d2_r,z2_r,x2_r,y2_r,age2_r,E_dot2_r,B2_r=[],[],[],[],[],[],[],[],[]
for i in range(len(P_dota_r)):
    if P_dota_r[i]!=0.0 and da_r[i]!=0.0 and E_dota_r[i]!=0 and da_r[i]<=25:
        P2_r+=[Pa_r[i]]
        P_dot2_r+=[P_dota_r[i]]
        d2_r+=[da_r[i]]
        z2_r+=[za_r[i]]
        x2_r+=[xa_r[i]]
        y2_r+=[ya_r[i]]
        age2_r+=[agea_r[i]]
        E_dot2_r+=[E_dota_r[i]]
        B2_r+=[Ba_r[i]]

for i in range(len(P_dota)):
    if P_dota[i]!=0.0 and da[i]!=0.0 and E_dota[i]!=0 and da[i]<=25:
        P2+=[Pa[i]]
        P_dot2+=[P_dota[i]]
        d2+=[da[i]]
        z2+=[za[i]]
        x2+=[xa[i]]
        y2+=[ya[i]]
        age2+=[agea[i]]
        E_dot2+=[E_dota[i]]
        B2+=[Ba[i]]

#Create file for optimisation program
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

#Log calculation for age, P and Pdot for the ATNF data
log_agea=[]
for i in range(len(age2)):
    log_agea+=[(np.log(agea[i]))/(np.log(10))]

log_Pa,log_P_dota=[],[]
for i in range(len(Pa)):
    log_Pa+=[np.log10(Pa[i])]
    log_P_dota+=[np.log10(P_dota[i])]

for i in range(len(P2)):
    log_P2+=[(np.log(P2[i]))/(np.log(10))]
    log_Pdot2+=[(np.log(P_dot2[i]))/(np.log(10))]

log_age_gamma_ATNF=[]
for i in range(len(age2_gall)):
    log_age_gamma_ATNF.append(np.log10(age2_gall[i]))

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

#with open("lat_dist_noISM.txt","w+")as f:
#    for i in range(len(latitude)):
#        f.write(f'{latitude[i]} {distance[i]} {x[i]} {y[i]}\n')

#Computation of characteristic age
charac_age,log_charac_age,P_old_charac=[],[],[]
for i in range(len(P)):
    charac_age.append((P[i]/(2*P_dot[i]))/(365*24*3600))
    log_charac_age.append(np.log10(charac_age[i]))

for i in range(len(P)):
    if charac_age[i]>=1e8:
        P_old_charac.append(np.log10(P[i]))

v0=[]
for i in range(len(vx0)):
    v0.append(np.sqrt(vx0[i]**2+vy0[i]**2+vz0[i]**2))

P_new,P_dot_new,age_new=[],[],[]
log_P_new,log_P_dot_new=[],[]
for i in range(len(P)):
    if age[i]<1e8*365*24*3600:
        P_new+=[P[i]]
        P_dot_new+=[P_dot[i]]
        age_new+=[age[i]]
        log_P_new+=[np.log10(P_new[i])]
        log_P_dot_new+=[np.log10(P_dot_new[i])]

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
P_gamma,P_dot_gamma,x_gamma,y_gamma,age_gamma,error_gamma,distance_gamma,cos_alpha_gamma=[],[],[],[],[],[],[],[]
P_radio_gamma,P_dot_radio_gamma,x_radio_gamma,y_radio_gamma,age_radio_gamma,error_radio_gamma,distance_radio_gamma=[],[],[],[],[],[],[]
log_age_gamma_all=[]
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
        cos_alpha_gamma+=[cos_alpha[i]]
    elif type_pulsar[i]==3:
        P_radio_gamma+=[P[i]]
        P_dot_radio_gamma+=[P_dot[i]]
        x_radio_gamma+=[x[i]]
        y_radio_gamma+=[y[i]]
        age_radio_gamma+=[age[i]]
        distance_radio_gamma+=[distance[i]]

P_gamma_or_rg,P_dot_gamma_or_rg=[],[]
P_gamma_or_rg=P_gamma+P_radio_gamma
P_dot_gamma_or_rg=P_dot_gamma+P_dot_radio_gamma

for i in range(len(age_gamma)):
    log_age_gamma_all.append(np.log10(age_gamma[i]/(365*24*3600)))

for i in range(len(age_radio_gamma)):
    log_age_gamma_all.append(np.log10(age_radio_gamma[i]/(365*24*3600)))

#Lists of the gamma-rays (gamma only or radio-gamma) detected pulsars above the usual threshold of Fermi/LAT
P_gamma_all,P_dot_gamma_all,x_gamma_all,y_gamma_all,age_gamma_all,error_gamma_all,distance_gamma_all,cos_alpha_gamma_all,latitude_gamma_all,longitude_gamma_all=[],[],[],[],[],[],[],[],[],[]
P_gamma_ab,P_dot_gamma_ab,x_gamma_ab,y_gamma_ab,age_gamma_ab,error_gamma_ab,distance_gamma_ab,cos_alpha_gamma_ab,latitude_gamma_ab,longitude_gamma_ab=[],[],[],[],[],[],[],[],[],[]
P_gamma_bel,P_dot_gamma_bel,x_gamma_bel,y_gamma_bel,age_gamma_bel,error_gamma_bel,distance_gamma_bel,cos_alpha_gamma_bel,latitude_gamma_bel,longitude_gamma_bel=[],[],[],[],[],[],[],[],[],[]
type_pulsar_g_or_rg=[]
for i in range(len(P)):
    if type_pulsar[i]==2 or type_pulsar[i]==3:
        P_gamma_all+=[P[i]]
        P_dot_gamma_all+=[P_dot[i]]
        x_gamma_all+=[x[i]]
        y_gamma_all+=[y[i]]
        age_gamma_all+=[age[i]]
        distance_gamma_all+=[distance[i]]
        cos_alpha_gamma_all+=[cos_alpha[i]]
        latitude_gamma_all+=[latitude[i]]
        longitude_gamma_all+=[longitude[i]]
        type_pulsar_g_or_rg+=[type_pulsar[i]]

for i in range(len(P_gamma_all)):
    if ((Fg_flux[i]<=16e-15 and type_pulsar_g_or_rg[i]==2) or (Fg_flux[i]<=4e-15 and type_pulsar_g_or_rg[i]==3)):
        P_gamma_ab+=[P_gamma_all[i]]
        P_dot_gamma_ab+=[P_dot_gamma_all[i]]
        x_gamma_ab+=[x_gamma_all[i]]
        y_gamma_ab+=[y_gamma_all[i]]
        age_gamma_ab+=[age_gamma_all[i]]
        distance_gamma_ab+=[distance_gamma_all[i]]
        cos_alpha_gamma_ab+=[cos_alpha_gamma_all[i]]
        latitude_gamma_ab+=[latitude_gamma_all[i]]
        longitude_gamma_ab+=[longitude_gamma_all[i]]
    elif ((Fg_flux[i]>16e-15 and type_pulsar_g_or_rg[i]==2) or (Fg_flux[i]>4e-15 and type_pulsar_g_or_rg[i]==3)):
        P_gamma_bel+=[P_gamma_all[i]]
        P_dot_gamma_bel+=[P_dot_gamma_all[i]]
        x_gamma_bel+=[x_gamma_all[i]]
        y_gamma_bel+=[y_gamma_all[i]]
        age_gamma_bel+=[age_gamma_all[i]]
        distance_gamma_bel+=[distance_gamma_all[i]]
        cos_alpha_gamma_bel+=[cos_alpha_gamma_all[i]]
        latitude_gamma_bel+=[latitude_gamma_all[i]]
        longitude_gamma_bel+=[longitude_gamma_all[i]]


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
P_line=[10**(np.log10(i)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line=[10**np.log10(const*(i**2)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line4=[Pdot_line[i]*10**(-0.55) for i in range(len(Pdot_line))]
Pdot_line5=[Pdot_line[i]*10**(1.15) for i in range(len(Pdot_line))]
Pdot_line2=[10**np.log10(const_min*(i**2)) for i in np.arange(1e-2,3e1,0.01)]
Pdot_line3=[10**np.log10(const_max*(i**2)) for i in np.arange(1e-2,3e1,0.01)]

#Plot the death line of the article of Chen & Ruderman (1993)
const_CR93=10**((43.8-16)/2)*16*(np.pi**3)*(R_NS**6)*(1+(np.sin(alpha_l)**2))/(Inertia*mu_0*(c_light**3))
P_dot_line_CR93=[10**np.log10(const_CR93*(i**(2))) for i in np.arange(1e-2,3e1,0.01)]

#Plot the line with the critical magnetic field 4.4e9 T, distinguishing magnetar from canonical pulsars
B_crit=4.4e9
B_linecrit=[10**np.log10(((B_crit/B0)**2)*(i**(-1))) for i in np.arange(1e-2,3e1,0.01)]

#Plot the Edot lines
P_dot_Edot1e22W=[(1e22*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e25W=[(1e25*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e28W=[(1e28*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]
P_dot_Edot1e31W=[(1e31*(P_line[i])**3)/(4*np.pi**2*Inertia) for i in range(len(P_line))]

#Plot the B lines
P_dot_B1e6=[(1e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B5e6=[(5e6/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e7=[(1e7/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e8=[(1e8/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B1e9=[(1e9/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]
P_dot_B4_4e9=[(4.4e9/3.2e15)**2*(1/P_line[i]) for i in range(len(P_line))]

#Display of the statistics on all pulsar
for i in range(len(Edot)):
    if Edot[i] > 1e31 and type_pulsar[i]==1:
        count_rad_E_bigbig+=1.0
    if Edot[i] > 1e31 and type_pulsar[i]==2:
        count_gam_E_bigbig+=1.0
    if Edot[i] > 1e31 and type_pulsar[i]==3:
        count_radgam_E_bigbig+=1.0
    if Edot[i] > 1e28 and type_pulsar[i]==1:
        count_rad_E_big+=1.0
    if Edot[i] > 1e28 and type_pulsar[i]==2:
        count_gam_E_big+=1.0
    if Edot[i] > 1e28 and type_pulsar[i]==3:
        count_radgam_E_big+=1.0
    if type_pulsar[i]==1:
        count_rad+=1.0
    if type_pulsar[i]==2:
        count_gam+=1.0
    if type_pulsar[i]==3:
        count_radgam+=1.0

print(f"Number of radio pulsars with Edot > 1e31 W : {count_rad_E_bigbig}\nNumber of gamma pulsars with Edot > 1e31 W : {count_gam_E_bigbig}")
print(f"Number of radio-gamma pulsars with Edot > 1e31 W : {count_radgam_E_bigbig}")
print(f"Number of radio pulsars with Edot > 1e28 W : {count_rad_E_big}")
print(f"Number of gamma pulsars with Edot > 1e28 W : {count_gam_E_big}")
print(f"Number of radio-gamma pulsars with Edot > 1e28 W : {count_radgam_E_big}")
print(f"Number of radio pulsars : {count_rad}")
print(f"Number of gamma pulsars : {count_gam}")
print(f"Number of radio-gamma pulsars : {count_radgam}")
count_tot=count_rad+count_gam+count_radgam
print(f"Number of total pulsars : {count_tot}")

#Make histograms to prepare for the 2D plot of comparison with observations
histSIM,xsim_edges,ysim_edges=np.histogram2d(log_P,log_Pdot,bins=(20,20))
histobs,xobs_edges,yobs_edges=np.histogram2d(log_Pa,log_P_dota,bins=(20,20))
histSIM=np.rot90(histSIM)
#histSIM=np.rot90(histSIM)
histobs=np.rot90(histobs)
#histobs=np.rot90(histobs)

hist_diff=((histSIM/len(log_P))-(histobs/len(log_Pa)))/(((histSIM/len(log_P))+(histobs/len(log_Pa)))**1.0)
chi2_tab=(histSIM-histobs)**2/histobs
chi2_tab=np.where(np.isfinite(chi2_tab), chi2_tab, 0)
chi2=np.nansum(chi2_tab)

#Make histograms to compute the CDF
counts,bin_edges=np.histogram(log_P,bins=20,density=False)
counts_2,bin_edges_2=np.histogram(log_Pdot,bins=20,density=False)
counts_3,bin_edges_3=np.histogram(log_P2,bins=20,density=False)
counts_4,bin_edges_4=np.histogram(log_Pdot2,bins=20,density=False)

cum_counts=np.cumsum(counts)
cum_counts_2=np.cumsum(counts_2)
cum_counts_3=np.cumsum(counts_3)
cum_counts_4=np.cumsum(counts_4)

total_count=cum_counts[-1]
total_count2=cum_counts_2[-1]
total_count3=cum_counts_3[-1]
total_count4=cum_counts_4[-1]

cdf1=cum_counts/total_count
cdf2=cum_counts_2/total_count2
cdf3=cum_counts_3/total_count3
cdf4=cum_counts_4/total_count4

#KS test 1D python (all the data)
KS_test=kstest(P_dot,P_dota)
test_stat_all1=KS_test.statistic
p_value_all1=KS_test.pvalue
print("----ALL THE DATA----\n")
print(f"d_value of Pdot KS test = {test_stat_all1}")
print(f"p_value of Pdot={p_value_all1}")

KS_test=kstest(P,Pa)
test_stat_all2=KS_test.statistic
p_value_all2=KS_test.pvalue
print(f"d_value of P KS test = {test_stat_all2}")
print(f"p_value of P={p_value_all2}")

#KS test 1D python (gamma-ray population)
KS_test=kstest(P_dot_gamma_or_rg,P_dota_gall)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL GAMMA PULSARS----\n")
print(f"d_value of Pdot KS test for all the gamma pulsars= {test_stat}")
print(f"p_value of Pdot for all the gamma pulsars={p_value}")

KS_test=kstest(P_gamma_or_rg,Pa_gall)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for all the gamma pulsars= {test_stat}")
print(f"p_value of P for all the gamma pulsars={p_value}")

#KS test 1D python (gamma-ray only population)
KS_test=kstest(P_dot_gamma,P_dota_gamma)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----GAMMA ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of Pdot for the gamma only pulsars={p_value}")

KS_test=kstest(P_gamma,Pa_gamma)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of P for the gamma only pulsars={p_value}")

#KS test 1D python (radio/gamma-ray population)
KS_test=kstest(P_dot_radio_gamma,P_dota_rg)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO/GAMMA PULSARS----\n")
print(f"d_value of Pdot KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of Pdot for the radio/gamma pulsars={p_value}")

KS_test=kstest(P_radio_gamma,Pa_rg)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of P for the radio/gamma pulsars={p_value}")

#KS test 1D python (radio population)
KS_test=kstest(P_dot_radio,P_dota_r)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the radio only pulsars= {test_stat}")
print(f"p_value of Pdot for the radio only pulsars={p_value}")

KS_test=kstest(P_radio,Pa_r)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio only pulsars= {test_stat}")
print(f"p_value of P for the radio only pulsars={p_value}")

#KS test without the outliers
KS_test=kstest(P_dot_gamma_or_rg,P_dota_gall_no5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL GAMMA PULSARS WITHOUT THE OUTLIERS----\n")
print(f"d_value of Pdot KS test for all the gamma pulsars= {test_stat}")
print(f"p_value of Pdot for all the gamma pulsars={p_value}")

KS_test=kstest(P_gamma_or_rg,Pa_gall_no5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for all the gamma pulsars= {test_stat}")
print(f"p_value of P for all the gamma pulsars={p_value}")

print(len(Pa_gall)-len(Pa_gall_no5))

#Plot the CDF
plt.plot(bin_edges[1:], cdf1, marker='o', linestyle='-',label='Simulation data')
plt.plot(bin_edges_3[1:], cdf3, marker='o', linestyle='-',label='ATNF data')
plt.xlabel('Log(P) (P in s)')
plt.ylabel('CDF of Log(P)')
plt.legend()
plt.savefig('CDF_log_P.png')
plt.close()

#Plot the CDF of log(P_dot)
plt.plot(bin_edges_2[1:], cdf2, marker='o', linestyle='-',label='Simulation data')
plt.plot(bin_edges_4[1:], cdf4, marker='o', linestyle='-',label='ATNF data')
plt.xlabel('Log(Pdot) (Pdot in s.s^-1)')
plt.ylabel('CDF of Log(Pdot)')
plt.legend()
plt.savefig('CDF_log_Pdot.png')
plt.close()

condition = [Pdot2 < Pdot3 for Pdot2, Pdot3 in zip(Pdot_line2,Pdot_line3)]

#P-Pdot plot only canonical population
plt.figure(300)
plt.scatter(Pa,P_dota,c='blue',marker='o',s=10,label='ATNF data')
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Mitra et al. 2019
plt.fill_between(P_line,Pdot_line4,Pdot_line5,where=condition,facecolor='green',alpha=0.4,label='Death Valley')

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

plt.legend(loc='lower left')
plt.xlim(1e-2,3e1)
plt.ylim(1e-19,1e-11)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio pulsars only")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_selected.png',dpi=300)
plt.close()

#P-Pdot plot all pulsars
plt.figure(1)
plt.scatter(P,P_dot,c='red',marker='o',s=5,label='Simulation data',zorder=2)
#plt.scatter(P2,P_dot2,c='blue',marker='o',s=5,label='ATNF data') #whole pop
plt.scatter(Pa,P_dota,c='blue',marker='o',s=5,label='ATNF data : canonical pulsars',zorder=1) #Only canonical pop
#plt.plot(P_line,Pdot_line4,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line5,c='green',linestyle='-',linewidth=2)
#plt.plot(P_line,Pdot_line3,c='brown')
#plt.plot(P_line,Pdot_line2,c='pink')
#plt.plot(P_death2,P_dot_death2,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Ruderman & Sutherland 1975 
plt.plot(P_line,Pdot_line,c='green',label='Death line',linestyle='-',linewidth=1) #Death line Mitra et al. 2019
#plt.plot(P_line,P_dot_line_CR93,c='green',label='Death line',linestyle='-',linewidth=2) #Death line Chen & Ruderman 1993
#plt.plot(P_line,B_linecrit,c='blue',label='Critical magnetic field line',linestyle='-',linewidth=2)
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
#plt.title("Spin period derivative - Spin period diagram")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.legend(loc='lower left',fontsize='x-small')
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
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
plt.savefig('2D P_Pdot_plot_sim.png',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(32)
#plt.figure(figsize=(6,6))
#plt.imshow((histSIM-histobs)/((histSIM+histobs)**0.5),extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.hist2d(log_P2,log_Pdot2,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
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
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
plt.savefig('2D P_Pdot_plot.png',dpi=300)
plt.close()

#P-Pdot plot radio pulsars only 
plt.figure(2)
plt.scatter(P_radio,P_dot_radio,c='red',marker='o',s=10,label='Simulation')
plt.scatter(Pa_r,P_dota_r,c='blue',marker='o',s=10,label='ATNF data')
plt.legend()
plt.xlim(1e-2,3e1)
plt.ylim(1e-18,1e-11)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio pulsars only")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_r.png')
plt.close()

#P-Pdot plot gamma pulsars only
plt.figure(3)
plt.scatter(P_gamma,P_dot_gamma,c='green',marker='o',s=5,label='Simulation')
plt.scatter(Pa_gamma,P_dota_gamma,c='blue',marker='o',s=5,label='ATNF data')
plt.xlim(1e-2,5)
plt.ylim(1e-17,1e-11)
plt.yscale('log')
plt.xscale('log')
plt.legend()
#plt.title("Spin period derivative - Spin period diagram for gamma pulsars only")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_gonly.png',dpi=300)
plt.close()

#P-Pdot plot radio gamma pulsars
plt.figure(4)
plt.scatter(P_radio_gamma,P_dot_radio_gamma,c='purple',marker='o',s=10,label='Simulation')
plt.scatter(Pa_rg,P_dota_rg,c='blue',marker='o',s=10,label='ATNF data')
plt.legend()
plt.xlim(1e-2,1e1)
plt.ylim(1e-20,1e-10)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio gamma pulsars")
plt.xlabel('Spin period s')
plt.ylabel('Spin period derivative s.s^-1')
plt.savefig('P_Pdot_plot_rg.png')
plt.close()

#P-Pdot plot all gamma pulsars
plt.figure(27)
plt.scatter(P_gamma_or_rg,P_dot_gamma_or_rg,c='green',marker='o',s=10,label='Simulation')
plt.scatter(Pa_gall,P_dota_gall,c='blue',marker='o',s=10,label='ATNF data')

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

plt.legend()
plt.xlim(1e-2,3e1)
plt.ylim(1e-19,1e-11)
plt.yscale('log')
plt.xscale('log')
#plt.title("Spin period derivative - Spin period diagram for radio gamma pulsars")
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_all_gamma.png',dpi=300)
plt.close()

#Positions plot
plt.figure(1000)
plt.scatter(x,y,s=2,c='red',alpha=0.4,label='Simulation with the ISM',zorder=2)
plt.scatter(x2,y2,s=2,c='blue',alpha=0.5,label='ATNF data') #Canonical pop
plt.scatter(y_dirson22,x_dirson22,s=2,c='green',alpha=0.7,label='Dirson et al.(2022)')
plt.scatter(x_noISM,y_noISM,s=2,c='black',alpha=0.4,label='Simulation without the ISM',zorder=0)
plt.scatter([0],[8.5],c='yellow',marker='o',s=20,label='The Sun',zorder=3) #position of the sun
plt.xlim(-20,20)
plt.ylim(-20,20)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.axis('equal')
#plt.title('Positions of the detected pulsars compared to the sun')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend(fontsize='small')
plt.savefig('Positions_detected_pulsars.png',dpi=300)
plt.close()

#Distance histogram
plt.figure(6)
#plt.hist([distance,d_selec,dist_dirson22],bins=20,range=(0,25),edgecolor='white',color=['red','blue','green'],label=['Simulation','ATNF data','Dirson et al. (2022)'])
plt.hist(distance,bins=20,range=(0,25),edgecolor='red',color='white',label='Simulation with the ISM',alpha=1,zorder=1,histtype='step')
#plt.hist(d2,bins=20,range=(0,25),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.hist(d2,bins=20,range=(0,25),edgecolor='blue',color='white',label='ATNF data',alpha=1,zorder=1,histtype='step') #Canonical pop
plt.hist(dist_dirson22,bins=20,range=(0,25),edgecolor='green',color='white',alpha=1,label='Dirson et al. (2022)',zorder=1,histtype='step')
plt.hist(dist_noISM,bins=20,range=(0,25),edgecolor='black',color='white',alpha=1,label='Simulation without the ISM',zorder=1,histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel('d (kpc)')
plt.ylabel('Frequency')
#plt.title('Distance to earth of the pulsars')
plt.savefig('histo_dist.png',dpi=300)
plt.close()

#Age histogram (charac age OBS VS real age SIM)
plt.figure(7)
plt.hist(log_age,bins=20,range=(1,11),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
#plt.hist(log_age2,bins=20,range=(1,11),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #Whole pop
plt.hist(log_agea,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
#plt.title('Histogram of the age of the detected pulsars')
plt.savefig('histo_age.png',dpi=300)
plt.close()

#Age histogram (charac age OBS VS charac age SIM)
plt.figure(100)
plt.hist(log_charac_age,bins=20,range=(1,11),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_agea,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
#plt.title('Histogram of the age of the detected pulsars')
plt.savefig('histo_age_both_charac.png',dpi=300)
plt.close()

#Longitude
plt.figure(27)
plt.hist(longitude,bins=20,range=(0,360),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
#plt.hist(latitude_selec,bins=20,range=(-60,60),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step') #Canonical pop
#plt.yscale('log')
plt.legend()
plt.xlabel('Longitude in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the latitude of the detected pulsars')
plt.savefig('histo_longitude.png',dpi=300)
plt.close()

#Latitude histogram
plt.figure(8)
plt.hist(latitude,bins=20,range=(-60,60),edgecolor='red',color='red',alpha=1,label='Simulation with the ISM',zorder=3,histtype='step')
plt.hist(latitude2,bins=20,range=(-60,60),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step') #Canonical pop
plt.hist(lat_dirson22,bins=20,range=(-60,60),edgecolor='green',color='green',alpha=1,label='Dirson et al.(2022)',zorder=1,histtype='step')
plt.hist(lat_noISM,bins=20,range=(-60,60),edgecolor='black',color='green',alpha=1,label='Simulation without the ISM',zorder=0,histtype='step')
plt.yscale('log')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
#plt.title('Histogram of the latitude of the detected pulsars')
plt.savefig('histo_latitude.png',dpi=300)
plt.close()

#Log(P) histogram
plt.figure(9)
plt.hist(log_P,bins=20,range=(-2,1.5),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_Pa,bins=20,range=(-2,1.5),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
#plt.hist(P_old_charac,bins=20,range=(-2,1.5),edgecolor='green',color='green',alpha=1,label='Old pulsars',histtype='step')
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period of the detected pulsars')
plt.savefig('histo_period.png',dpi=300)
plt.close()

#Log(P_dot) histogram
plt.figure(10)
plt.hist(log_Pdot,bins=20,range=(-20,-10),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_P_dota,bins=20,range=(-20,-10),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel(r'Log ($\dot{P}$) ($\ \dot{P}$ in $s.s^{-1})$')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period derivative of the detected pulsars')
plt.savefig('histo_pdot.png',dpi=300)
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(11)
plt.hist(cos_alpha,bins=20,range=(-1,1),edgecolor='black',color='red',alpha=0.5,label=r'cos($\alpha$)')
plt.hist(cos_alpha0,bins=20,range=(-1,1),edgecolor='black',color='blue',alpha=0.5,label=r'cos($\alpha_0$)')
plt.legend()
plt.xlabel(r'cos($\alpha$) and cos($\alpha_0$)')
plt.ylabel('Frequency')
#plt.title('Histogram of the inclination angle of the detected pulsars')
plt.savefig('histo_cosalpha.png')
plt.close()

#Period of old pulsars : histogram
plt.figure(12)
plt.hist(P_old,bins=20,range=(-2,1.5),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('Frequency')
#plt.title('Histogram of the rotation period of the detected pulsars with ages betwen 1e7.7 years and 1e8.7')
plt.savefig('histo_period_old.png')
plt.close()

#Relative error between B(P,Pdot) and Bf the decaying magnetic field
plt.figure(13)
plt.hist(err_rel_B,bins=20,range=(0,0.01),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel(r'Relative error between the magnetic field computed with decay or with $P$ and $\dot{P}$')
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
plt.hist(B_ppdot,bins=20,range=(1e6,1.5e8),edgecolor='black',color='blue',alpha=0.5,label=r'Magnetic field computed with $P$ and $\dot{P}$')
plt.legend()
plt.xlabel('Magnetic field (B in Tesla)')
plt.ylabel('Frequency')
plt.savefig('histo_Bfield.png')
plt.close()

#Spin-velocity angles histogram
plt.figure(16)
plt.hist(PA,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Frequency')
plt.yscale('log')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle.png',dpi=300)
plt.close()

#Spin-velocity angle (old pulsars) histogram
plt.figure(17)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees for pulsars older than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_old.png')
plt.close()

#Spin-velocity angle (young pulsars) histogram
plt.figure(18)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees for pulsars younger than 10 Myr')
plt.ylabel('Frequency')
#plt.title('Histogram of the angle between the velocity vector and the rotation axis of the detected pulsars')
plt.savefig('histo_spinvelangle_young.png')
plt.close()

#Spin-velocity angle (young pulsars and old) histogram
plt.figure(19)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=1,label='Age < 10Myr',zorder=2,histtype='step')
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Age > 10Myr',histtype='step')
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

#cos(alpha0) and cos(alpha) gamma pulsars histogram
plt.figure(22)
plt.hist(cos_alpha_gamma,bins=20,range=(-1,1),edgecolor='black',color='green',alpha=0.5,label=r'cos($\alpha$)')
plt.legend()
plt.xlabel(r'cos($\alpha$)')
plt.ylabel('Frequency')
plt.savefig('histo_cosalpha_gamma_ray_pulsars.png')
plt.close()

#Age histogram every gamma pulsars 
plt.figure(23)
plt.hist(log_age_gamma_all,bins=20,range=(1,11),edgecolor='green',color='green',alpha=1,label='Simulation',histtype='step')
plt.hist(log_age2,bins=20,range=(1,11),edgecolor='black',color='blue',alpha=0.5,label='ATNF data') #canonical pop
plt.hist(log_age_gamma_ATNF,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #gamma pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Frequency')
plt.savefig('histo_age_gamma.png',dpi=300)
plt.close()

#P-Pdot every gamma pulsars with threshold below the Fermi/LAT sensitivity
plt.figure(24)
plt.scatter(P_gamma_ab,P_dot_gamma_ab,c='green',marker='o',s=5,label='Simulation')
plt.scatter(P2_gall,P_dot2_gall,c='blue',marker='o',s=5,label='ATNF data')
plt.xlim(1e-2,3e1)
plt.ylim(1e-18,1e-11)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_g_above_Fermi_LAT.png',dpi=300)
plt.close()

#Latitude gamma pulsars below threshold sensitivity of Fermi/LAT
plt.figure(25)
plt.hist(latitude_gamma_ab,bins=20,range=(-60,60),edgecolor='green',color='green',alpha=1,label='Simulation',zorder=3,histtype='step')
plt.yscale('log')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Frequency')
plt.savefig('histo_latitude_above_Fermi_LAT.png',dpi=300)
plt.close()

print(len(x_gamma_ab))
print(len(x_gamma_bel))

#X-Y positions gamma pulsars below threshold sensitivity of Fermi/LAT
plt.figure(26)
plt.scatter(x_gamma_ab,y_gamma_ab,s=2,c='green',alpha=0.2,label='Detection with a sensitivity lower than a factor of 10 compared to\nthe Fermi/LAT instrument',zorder=2)#label=r'Gamma pulsars with $F_{min} < 16\times10^{-15} W.m^{-2}$ or $F_{min} < 4\times10^{-15} W.m^{-2}$',zorder=2)
plt.scatter(x_gamma_bel,y_gamma_bel,s=2,c='red',alpha=0.9,label='Detection according to the Fermi/LAT instrument sensitivity',zorder=1)#label=r'Gamma pulsars with $F_{min} > 16\times10^{-15} W.m^{-2}$ or $F_{min} > 4\times10^{-15} W.m^{-2}$',zorder=1)
plt.scatter([0],[8.5],c='yellow',marker='o',s=20,label='The Sun',zorder=3) #position of the sun
plt.xlim(-15,15)
plt.ylim(-15,15)
#plt.gca().set_aspect('equal', adjustable='box')
#plt.axis('equal')
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend(fontsize='small')
plt.savefig('Positions_detected_pulsars_above_Fermi_LAT.png',dpi=300)
plt.close()

for i in range(len(age)):
    age[i]=age[i]/(365*24*3600) #Get the age in yr

#Plot spin-vel angle=f(nb_orbit)
plt.figure(40)
scatter=plt.scatter(PA,nb_orbit,c=age,cmap='plasma',marker='o',s=5,label='Simulation data')#,norm=LogNorm(vmin=1e2,vmax=5e7))
#plt.ylim(0,10)
plt.xlim(0,180)
plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Number of orbits in the Galaxy')
plt.legend()
cbar=plt.colorbar(scatter)
cbar.set_label('Age (yr)')
plt.savefig('spinvel_orbit_plot.png',dpi=300)
plt.close()

#Plot spin-vel angle=f(v_init)
plt.figure(41)
scatter=plt.scatter(PA,nb_orbit,c=v0,cmap='plasma',marker='o',s=5,label='Simulation data')
plt.ylim(0,10)
plt.xlim(0,180)
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Number of orbits in the Galaxy')
cbar=plt.colorbar(scatter)
cbar.set_label('Initial velocity (km/s)')
plt.legend()
plt.savefig('spinvel_init_vel_plot.png',dpi=300)
plt.close()

Fg_flux_log=[]
for i in range(len(Fg_flux)):
    Fg_flux_log.append(np.log10(Fg_flux[i]))

flux_3PC_cano_log=[]
for i in range(len(flux_3PC_cano)):
    flux_3PC_cano_log.append(np.log10(flux_3PC_cano[i]))

#Comparison flux 3PC and simulation (gamma-ray pulsars)
plt.figure(42)
plt.hist(Fg_flux_log,bins=20,range=(-15,-11),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
plt.hist(flux_3PC_cano_log,bins=20,range=(-15,-11),edgecolor='blue',color='blue',alpha=1,label='3PC data',zorder=2,histtype='step') #Canonical pop
plt.legend()
plt.yscale('log')
plt.xlabel(r'Gamma-ray flux in log space (W.m$^{-2}$)')
plt.ylabel('Frequency')
plt.savefig('histo_fgamma.png',dpi=300)
plt.close()

weights_obs=np.ones_like(g_peak_sep_obs) / len(g_peak_sep_obs)
weights_sim=np.ones_like(g_peak_sep_sim) / len(g_peak_sep_sim)

#Comparison gamma-ray peak separation 3PC and simulation
plt.figure(43)
plt.hist(g_peak_sep_sim,bins=10,weights=weights_sim,range=(0,0.5),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
plt.hist(g_peak_sep_obs,bins=10,weights=weights_obs,range=(0,0.5),edgecolor='blue',color='blue',alpha=1,label='3PC data',zorder=2,histtype='step') #Canonical pop
plt.legend()
#plt.yscale('log')
plt.xlabel(r'Gamma-ray peak separation')
plt.ylabel('p.d.f value')
plt.savefig('histo_gpeaksep.png',dpi=300)
plt.close()

#KS test without the outliers
KS_test=kstest(Fg_flux_log,flux_3PC_cano_log)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL GAMMA PULSARS FLUX----\n")
print(f"d_value of Fg KS test for all the gamma pulsars= {test_stat}")
print(f"p_value of Fg for all the gamma pulsars={p_value}")

#"Optimisation save in file" 
#line=[f"{count_tot} {count_rad} {count_gam} {count_radgam} {test_stat_all2} {p_value_all2} {test_stat_all1} {p_value_all1}\n"]
#with open("data_opt_04sigp.txt","a") as f:
#    f.writelines(line)
