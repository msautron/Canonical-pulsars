import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d
from scipy.stats import kstest,ks_2samp,linregress
from scipy.stats import mannwhitneyu
from matplotlib.colors import LogNorm
from astropy.table import Table
import pandas as pd

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out 
Pa_gamma,P_dota_gamma,da_gamma,za_gamma,xa_gamma,ya_gamma,agea_gamma,E_dota_gamma=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue of the gamma pulsars only
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used
P3,P_dot3,d3,z3,x3,y3,age3,E_dot3,latitude3,longitude3=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio)
P4,P_dot4,d4,z4,x4,y4,age4,E_dot4,latitude4,longitude4=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (gamma)
P5,P_dot5,d5,z5,x5,y5,age5,E_dot5,latitude5,longitude5=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio-gamma)
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
test_l=[]
B4,B5,B3=[],[],[]
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
reg_5=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
reg_6=re.compile("-?\d+[.]\d+")
reg_survey=re.compile(r"\b(?!NULL\b)[\w,]{3,}\b")

with open("fast_fermi_pmps.txt","r") as f:
    data2=re.findall(reg_3,f.read())

with open("fast_fermi_pmps.txt","r") as f:
    data_survey=re.findall(reg_survey,f.read())

with open("w10_gl_gb.txt","r") as f:
    data_gl_gb_w10_ATNF=re.findall(reg_6,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data=re.findall(reg_1,f.read())

with open("P_Pdot_positions.txt","r") as f:
    data_type=re.findall(reg_2,f.read())

with open("Fg_flux.txt","r") as f:
    fg_data=re.findall(reg_4,f.read())

with open("gamma_peak_sep.txt","r") as f:
    g_peak_sep_data=re.findall(reg_4,f.read())

with open("nb_orbit.txt","r") as f:
    nb_orbit_data=re.findall(reg_4,f.read())

with open("wint.txt","r") as f:
    data_wint=re.findall(reg_5,f.read())

with open("xi_rho_data.txt","r") as f:
    xi_rho_data=re.findall(reg_5,f.read())

with open("wr.txt","r") as f:
    data_wr=re.findall(reg_5,f.read())

df=pd.read_excel('3PC_Catalog_20230803.xls')
data_3PC=Table.from_pandas(df)

#Get the flux from the 3PC catalog for the canonical pulsars
flux_3PC_cano,g_peak_sep_obs=[],[]
count_msp_g=0
for i in range(len(data_3PC['G100'])):
    if data_3PC['B_S'][i]!='*':
        if float(data_3PC['B_S'][i])<=6e10:
            count_msp_g+=1
        if float(data_3PC['B_S'][i])<=4.4e13 and float(data_3PC['B_S'][i])>=6e10 and data_3PC['G100'][i]!='*':
            flux_3PC_cano.append(float(data_3PC['G100'][i])*1e-3)
        if float(data_3PC['B_S'][i])<=4.4e13 and float(data_3PC['B_S'][i])>=6e10 and data_3PC['PKSEP'][i]!='*':
            if data_3PC['PKSEP'][i] <= 0.5: 
                g_peak_sep_obs.append(float(data_3PC['PKSEP'][i]))
            else:
                g_peak_sep_obs.append(1.0-float(data_3PC['PKSEP'][i]))

print(f'Number of canonical gamma pulsars : {len(g_peak_sep_obs)}')
print(f'Number of gamma MSP in 3PC : {count_msp_g}')

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

#Get the data from the ATNF catalog for longitude and latitude
latitude_ATNF,longitude_ATNF,w10_atnf=[],[],[]
for i in range(int(len(data_gl_gb_w10_ATNF)/3)):
    longitude_ATNF.append(float(data_gl_gb_w10_ATNF[3*i]))
    latitude_ATNF.append(float(data_gl_gb_w10_ATNF[3*i+1]))
    w10_atnf.append(float(data_gl_gb_w10_ATNF[3*i+2]))

#Get the data from the wr file for each MSPs
wr_cano=[]
for i in range(int(len(data_wr))):
        wr_cano.append(float(data_wr[i]))

wint_cano=[]
for i in range(int(len(data_wint))):
    wint_cano.append(float(data_wint[i]))

#Get rho et xi from simulation
rho_sim,xi_sim=[],[]
for i in range(int(len(xi_rho_data)/2)):
    xi_sim.append(float(xi_rho_data[2*i])*180/np.pi)
    rho_sim.append(float(xi_rho_data[2*i+1])*180/np.pi)

#Get data from 3PC and filter the MSPs
df=pd.read_excel('3PC_Catalog_20230803.xls')
data_3PC=Table.from_pandas(df)
data_3PC['B_S']=pd.to_numeric(data_3PC['B_S'],errors='coerce')
mask=(data_3PC['B_S'] > 6e10) & (data_3PC['B_S'] < 4.4e13)
data_3PC_filtered=data_3PC[mask]

size_3PC_filtered=len(data_3PC_filtered['P0'])

print(f'Number of canonical pulsars in 3PC {size_3PC_filtered}')

#name_to_exclude=['J0023+0923','J0340+4130','J0605+3757','J0653+4706','J1125-6014','J1142+0119','J1301+0833','J1312+0051','J1455-3330','J1544+4937','J1630+3734','J1730-2304','J1741+1351','J1745+1017','J1824-2452A','J1843-1113','J1921+0137','J1921+1929','J1946+3417','J1959+2048','J2017+0603','J2042+0246','J2215+5135','J2234+0944']
#mask_exclusion = ~np.isin(data_3PC_filtered['PSRJ'], name_to_exclude)
#data_3PC_filtered_complete = data_3PC_filtered[mask_exclusion]

#Handle the type of pulsar
type_pulsar_obs=[]
for i in range(len(data_survey)):
    if ("fast" in data_survey[i] or "pks" in data_survey[i]) and not "Fermi" in data_survey[i]:
        type_pulsar_obs.append(1) #Pulsar radio
    elif not ("fast" in data_survey[i] and not "pks" in data_survey[i]) and "Fermi" in data_survey[i]:
        type_pulsar_obs.append(2) #Pulsar gamma
    elif ("fast" in data_survey[i] or "pks" in data_survey[i]) and "Fermi" in data_survey[i]:
        type_pulsar_obs.append(3) #Pulsar rg

#Get the data of the ATNF catalogue
for i in range(int(len(data2)/9)):
    Pa+=[float(data2[i*9])] #Period of the rotation of the pulsar in seconds
    P_dota+=[float(data2[i*9+1])] #Period derivative of the rotation of the pulsar no units
    da+=[float(data2[i*9+2])] #Distance to us in kpc
    za+=[float(data2[i*9+3])] #Z position in the galactocentric frame in kpc
    xa+=[float(data2[i*9+4])] #X position in the galactocentric frame in kpc
    ya+=[float(data2[i*9+5])] #age of the pulsar in years
    agea+=[float(data2[i*9+6])] #age of the pulsars in seconds
    Ba+=[float(data2[i*9+7])*1e-4] #surface magnetic field in T
    E_dota+=[float(data2[i*9+8])] #Spin down power of the pulsar in ergs/s

#Data we are really using, getting rid of lines where we miss info
type_pulsar_obs2,w10_atnf2=[],[]
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
        w10_atnf2.append(w10_atnf[i])
        type_pulsar_obs2.append(type_pulsar_obs[i])

#Computation of latitude and longitude for the ATNF data
for i in range(len(z2)):
    if d2[i]<1e-15:
        lat=0
        latitude2+=[lat]
    else:
        lat=np.arcsin(z2[i]/d2[i])*RAD
        latitude2+=[lat]
    r=(x2[i]**2+y2[i]**2)**0.5
    if x2[i]>=0:
        longi=np.arccos(-y2[i]/r)*RAD
        longitude2+=[longi]
    else:
        longi=(np.arccos(-y2[i]/r)+np.pi)*RAD
        longitude2+=[longi]

w10_r,w10_rg=[],[]
for i in range(len(type_pulsar_obs2)):
    if type_pulsar_obs2[i]==1:
        P3+=[P2[i]]
        P_dot3+=[P_dot2[i]]
        d3+=[d2[i]]
        z3+=[z2[i]]
        x3+=[x2[i]]
        y3+=[y2[i]]
        age3+=[age2[i]]
        E_dot3+=[E_dot2[i]]
        B3+=[B2[i]]
        latitude3+=[latitude2[i]]
        longitude3+=[longitude2[i]]
        w10_r+=[(float(w10_atnf2[i])*(1e-3)*(360/P2[i]))%360]
    elif type_pulsar_obs2[i]==2:
        P4+=[P2[i]]
        P_dot4+=[P_dot2[i]]
        d4+=[d2[i]]
        z4+=[z2[i]]
        x4+=[x2[i]]
        y4+=[y2[i]]
        age4+=[age2[i]]
        E_dot4+=[E_dot2[i]]
        B4+=[B2[i]]
        latitude4+=[latitude2[i]]
        longitude4+=[longitude2[i]]
    elif type_pulsar_obs2[i]==3:
        P5+=[P2[i]]
        P_dot5+=[P_dot2[i]]
        d5+=[d2[i]]
        z5+=[z2[i]]
        x5+=[x2[i]]
        y5+=[y2[i]]
        age5+=[age2[i]]
        E_dot5+=[E_dot2[i]]
        B5+=[B2[i]]
        latitude5+=[latitude2[i]]
        longitude5+=[longitude2[i]]
        w10_rg+=[(float(w10_atnf2[i])*(1e-3)*(360/P2[i]))%360]

print(f'Number of observed radio pulsars (FAST GPPS + PMPS + Arecibo): {len(P3)}')
print(f'Number of observed gamma pulsars (Fermi): {len(P4)}')
print(f'Number of observed radio+gamma pulsars: {len(P5)}')

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

#Refolding wr
for i in range(len(wr_cano)):
    if wr_cano[i]>360:
        wr_cano[i]=wr_cano[i]-int(wr_cano[i]/360.0)*360

#Lists depending on the pulsar emission type
P_radio,P_dot_radio,x_radio,y_radio,age_radio,error_radio,distance_radio=[],[],[],[],[],[],[]
P_gamma,P_dot_gamma,x_gamma,y_gamma,age_gamma,error_gamma,distance_gamma=[],[],[],[],[],[],[]
P_radio_gamma,P_dot_radio_gamma,x_radio_gamma,y_radio_gamma,age_radio_gamma,error_radio_gamma,distance_radio_gamma=[],[],[],[],[],[],[]
wr_r_or_rg,P_r_or_rg=[],[]
xi_r,xi_g,xi_rg=[],[],[]
alpha_r,alpha_g,alpha_rg=[],[],[]
Bf_r,Bf_g,Bf_rg=[],[],[]
z_r,z_g,z_rg=[],[],[]
Fr_gamma=[]
Edot_g=0
w_geometry_r_or_rg=[]
rho_r,rho_g,rho_rg=[],[],[]
xi_r_or_rg,rho_r_or_rg,alpha_r_or_rg=[],[],[]
for i in range(len(P)):
    if type_pulsar[i]==1:
        P_radio+=[P[i]]
        P_dot_radio+=[P_dot[i]]
        x_radio+=[x[i]]
        y_radio+=[y[i]]
        age_radio+=[age[i]]
        distance_radio+=[distance[i]]
        wr_r_or_rg.append(wr_cano[i])
        w_geometry_r_or_rg.append(wint_cano[i])
        P_r_or_rg.append(P[i])
        xi_r.append(xi_sim[i])
        xi_r_or_rg.append(xi_sim[i])
        alpha_r.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))
        alpha_r_or_rg.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))
        Bf_r.append(Bf[i])
        z_r.append(z[i])
        rho_r.append(rho_sim[i])
        rho_r_or_rg.append(rho_sim[i])
    elif type_pulsar[i]==2:
        P_gamma+=[P[i]]
        P_dot_gamma+=[P_dot[i]]
        x_gamma+=[x[i]]
        y_gamma+=[y[i]]
        age_gamma+=[age[i]]
        distance_gamma+=[distance[i]]
        xi_g.append(xi_sim[i])
        alpha_g.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))
        Bf_g.append(Bf[i])
        z_g.append(z[i])
        Edot_g=4*np.pi**2*1e38*(P_dot[i])*(P[i])**(-3)
        Fj=np.random.normal(loc=0.0,scale=0.2)
        Fr_gamma_value=9*1e3*(distance[i])**(-2)*(Edot_g/1e29)**(0.25)*10**Fj
        Fr_gamma.append(Fr_gamma_value)
        rho_g.append(rho_sim[i])
    elif type_pulsar[i]==3:
        P_radio_gamma+=[P[i]]
        P_dot_radio_gamma+=[P_dot[i]]
        x_radio_gamma+=[x[i]]
        y_radio_gamma+=[y[i]]
        age_radio_gamma+=[age[i]]
        distance_radio_gamma+=[distance[i]]
        wr_r_or_rg.append(wr_cano[i])
        w_geometry_r_or_rg.append(wint_cano[i])
        P_r_or_rg.append(P[i])
        xi_rg.append(xi_sim[i])
        xi_r_or_rg.append(xi_sim[i])
        alpha_rg.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))
        alpha_r_or_rg.append(min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi))
        Bf_rg.append(Bf[i])
        z_rg.append(z[i])
        rho_rg.append(rho_sim[i])
        rho_r_or_rg.append(rho_sim[i])

#Checking condition of radio observation
Pcheck,Pdotcheck=[],[]
xicheck,rhocheck,alphacheck,alpha_xi=[],[],[],[]
for i in range(len(P_r_or_rg)):
    if (wr_r_or_rg[i]<1):
        Pcheck.append(P_r_or_rg[i])
        xicheck.append(xi_r_or_rg[i])
        rhocheck.append(rho_r_or_rg[i])
        alphacheck.append(alpha_r_or_rg[i])
        alpha_xi.append(alpha_r_or_rg[i]+xi_r_or_rg[i])

#Compute the number of pulsars in the Galactic center
dist_GC=[]
nb_in_GC=0
L_gamma_all_inGC=0
Edot_calc=0
for i in range(len(x)):
    dist_calc=np.sqrt((x[i])**2+(y[i])**2+(z[i])**2)
    if (dist_calc <= 0.7):
        Edot_calc=4*np.pi**2*Inertia*P[i]**(-3)*P_dot[i]
        L_gamma_all_inGC+=10**(26.15)*(Bf[i]/1e8)**(0.11)*(Edot_calc/1e26)**(0.51)
        nb_in_GC+=1

print(f'Number of pulsars in the GC (detected or not) : {nb_in_GC}\nLuminosity of pulsars in the GC : {L_gamma_all_inGC}\n')

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

#KS test 1D python (gamma-ray only population)
KS_test=kstest(P_dot_gamma,data_3PC_filtered['P1'])
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----GAMMA ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of Pdot for the gamma only pulsars={p_value}")

KS_test=kstest(P_gamma,data_3PC_filtered['P0'])
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the gamma only pulsars= {test_stat}")
print(f"p_value of P for the gamma only pulsars={p_value}")

#KS test 1D python (radio/gamma-ray population)
KS_test=kstest(P_dot_radio_gamma,P_dot5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO/GAMMA PULSARS----\n")
print(f"d_value of Pdot KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of Pdot for the radio/gamma pulsars={p_value}")

KS_test=kstest(P_radio_gamma,P5)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio/gamma pulsars= {test_stat}")
print(f"p_value of P for the radio/gamma pulsars={p_value}")

#KS test 1D python (radio population)
KS_test=kstest(P_dot_radio,P_dot3)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print("----ALL RADIO ONLY PULSARS----\n")
print(f"d_value of Pdot KS test for the radio only pulsars= {test_stat}")
print(f"p_value of Pdot for the radio only pulsars={p_value}")

KS_test=kstest(P_radio,P3)
test_stat=KS_test.statistic
p_value=KS_test.pvalue
print(f"d_value of P KS test for the radio only pulsars= {test_stat}")
print(f"p_value of P for the radio only pulsars={p_value}")

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
plt.figure(0)
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
plt.savefig('P_Pdot_plot_selected.pdf',dpi=300)
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
plt.savefig('P_Pdot_plot.pdf',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(2)
#plt.figure(figsize=(6,6))
#plt.imshow((histSIM-histobs)/((histSIM+histobs)**0.5),extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(xsim_edges[:-1],ysim_edges[:-1],bins=(xsim_edges,ysim_edges),cmap='RdBu',weights=hist_diff.flatten())
plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
plt.savefig('2D P_Pdot_plot_sim.pdf',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(3)
#plt.figure(figsize=(6,6))
#plt.imshow((histSIM-histobs)/((histSIM+histobs)**0.5),extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.hist2d(log_P2,log_Pdot2,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
plt.savefig('2D P_Pdot_plot_obs.pdf',dpi=300)
plt.close()

#2D histogram P-Pdot
plt.figure(4)
#plt.figure(figsize=(6,6))
plt.imshow(hist_diff,extent=[-2, 1.5, -18.5, -11], cmap='RdBu',aspect='auto')
#plt.hist2d(xsim_edges[:-1],ysim_edges[:-1],bins=(xsim_edges,ysim_edges),cmap='RdBu',weights=hist_diff.flatten())
#plt.hist2d(log_P,log_Pdot,bins=(30,30),range=((-2,1.5),(-20,-10)),cmap='RdBu')
plt.colorbar()
#plt.gca().set_aspect('equal')
plt.xlabel(r'Log $P$ ($P$ in s)')
plt.ylabel(r'Log $\dot{P}$ ($\ \dot{P}$ in $s.s^{-1})$')
plt.savefig('2D P_Pdot_plot.pdf',dpi=300)
plt.close()

#P-Pdot plot radio pulsars only 
plt.figure(5)
plt.scatter(P_radio,P_dot_radio,c='red',marker='o',s=10,label='Simulation')
plt.scatter(P3,P_dot3,c='blue',marker='o',s=10,label='ATNF data')
#Plot the Edot lines
plt.plot(P_line,P_dot_Edot1e23W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[6]*1.0,P_dot_Edot1e23W[6]*1.5,r'$\dot{E} = 10^{23}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[15]*1.0,P_dot_Edot1e28W[15]*1.5,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
plt.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
plt.text(P_line[1]*0.68,P_dot_Edot1e25W[1]*0.5,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[60]*0.6,P_dot_B5e6[50]*0.64,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
#plt.plot(P_line,P_dot_B1e3,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[100]*0.6,P_dot_B1e3[90]*0.95,r'$B= 10^{3}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e4,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[800]*0.6,P_dot_B1e4[750]*0.95,r'$B= 10^{4}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e5,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e5[750]*0.9,r'$B= 10^{5}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
plt.text(P_line[900]*0.6,P_dot_B1e6[750]*0.9,r'$B= 10^{6}$T',fontsize=7,color='black',rotation=-17,zorder=3)
plt.legend()
plt.xlim(1e-2,3e1)
plt.ylim(1e-18,1e-11)
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_r.pdf')
plt.close()

#P-Pdot plot gamma pulsars only
plt.figure(6)
plt.scatter(P_gamma,P_dot_gamma,c='green',marker='o',s=5,label='Simulation')
plt.scatter(data_3PC_filtered['P0'],data_3PC_filtered['P1'],c='blue',marker='o',s=5,label='3PC data')
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
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_gonly.pdf',dpi=300)
plt.close()

#P-Pdot plot radio gamma pulsars
plt.figure(7)
plt.scatter(P_radio_gamma,P_dot_radio_gamma,c='purple',marker='o',s=10,label='Simulation')
plt.scatter(P5,P_dot5,c='blue',marker='o',s=10,label='ATNF data')
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
plt.legend()
plt.xlabel(r'$P$ (s)')
plt.ylabel(r'$\dot{P} \ (s.s^{-1})$')
plt.savefig('P_Pdot_plot_rg.pdf')
plt.close()

#Positions plot
plt.figure(8)
plt.scatter(x,y,s=2,c='red',alpha=0.4,label='Simulation with the ISM',zorder=2)
plt.scatter(x2,y2,s=2,c='blue',alpha=0.5,label='ATNF data') #Canonical pop
plt.scatter([0],[8.5],c='yellow',marker='o',s=20,label='The Sun',zorder=3) #position of the sun
plt.xlim(-20,20)
plt.ylim(-20,20)
plt.xlabel('x (kpc)')
plt.ylabel('y (kpc)')
plt.legend(fontsize='small')
plt.savefig('Positions_detected_pulsars.pdf',dpi=300)
plt.close()

#Distance histogram
plt.figure(9)
plt.hist(distance,bins=20,range=(0,25),edgecolor='red',color='white',label='Simulation with the ISM',alpha=1,zorder=1,histtype='step')
plt.hist(d2,bins=20,range=(0,25),edgecolor='blue',color='white',label='ATNF data',alpha=1,zorder=1,histtype='step') #Canonical pop
plt.legend()
plt.yscale('log')
plt.xlabel('d (kpc)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_dist.pdf',dpi=300)
plt.close()

#Age histogram (charac age OBS VS real age SIM)
plt.figure(10)
plt.hist(log_age,bins=20,range=(1,11),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_agea,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_age.pdf',dpi=300)
plt.close()

#Age histogram (charac age OBS VS charac age SIM)
plt.figure(11)
plt.hist(log_charac_age,bins=20,range=(1,11),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_agea,bins=20,range=(1,11),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel('Log(age) (age in yr)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_age_both_charac.pdf',dpi=300)
plt.close()

#Longitude
plt.figure(12)
plt.hist(longitude,bins=20,range=(0,360),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
#plt.hist(latitude_selec,bins=20,range=(-60,60),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step') #Canonical pop
#plt.yscale('log')
plt.legend()
plt.xlabel('Longitude in degrees')
plt.ylabel('Number of pulsars')
plt.savefig('histo_longitude.pdf',dpi=300)
plt.close()

#Latitude histogram
plt.figure(13)
plt.hist(latitude,bins=20,range=(-60,60),edgecolor='red',color='red',alpha=1,label='Simulation with the ISM',zorder=3,histtype='step')
plt.hist(latitude2,bins=20,range=(-60,60),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step') #Canonical pop
plt.yscale('log')
plt.legend()
plt.xlabel('Latitude in degrees')
plt.ylabel('Number of pulsars')
plt.savefig('histo_latitude.pdf',dpi=300)
plt.close()

#Log(P) histogram
plt.figure(14)
plt.hist(log_P,bins=20,range=(-2,1.5),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_Pa,bins=20,range=(-2,1.5),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_period.pdf',dpi=300)
plt.close()

#Log(P_dot) histogram
plt.figure(15)
plt.hist(log_Pdot,bins=20,range=(-20,-10),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.hist(log_P_dota,bins=20,range=(-20,-10),edgecolor='blue',color='blue',alpha=1,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel(r'Log ($\dot{P}$) ($\ \dot{P}$ in $s.s^{-1})$')
plt.ylabel('Number of pulsars')
plt.savefig('histo_pdot.pdf',dpi=300)
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(16)
plt.hist(cos_alpha,bins=20,range=(-1,1),edgecolor='black',color='red',alpha=0.5,label=r'cos($\alpha$)')
plt.hist(cos_alpha0,bins=20,range=(-1,1),edgecolor='black',color='blue',alpha=0.5,label=r'cos($\alpha_0$)')
plt.legend()
plt.xlabel(r'cos($\alpha$) and cos($\alpha_0$)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_cosalpha.pdf')
plt.close()

#Period of old pulsars : histogram
plt.figure(17)
plt.hist(P_old,bins=20,range=(-2,1.5),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_period_old.pdf')
plt.close()

#Relative error between B(P,Pdot) and Bf the decaying magnetic field
plt.figure(18)
plt.hist(err_rel_B,bins=20,range=(0,0.01),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel(r'Relative error between the magnetic field computed with decay or with $P$ and $\dot{P}$')
plt.ylabel('Number of pulsars')
plt.savefig('histo_err_B.pdf')
plt.close()

#Velocities
plt.figure(19)
plt.hist(v_norm,bins=20,range=(0,1000),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Velocities of the pulsars in km/s')
plt.ylabel('Number of pulsars')
plt.savefig('histo_velocities.pdf')
plt.close()

#Comparison between Bfield computed with the decay and with P and Pdot
plt.figure(20)
plt.hist(Bf,bins=20,range=(1e6,1.5e8),edgecolor='black',color='red',alpha=0.5,label='Decaying magnetic field')
plt.hist(B_ppdot,bins=20,range=(1e6,1.5e8),edgecolor='black',color='blue',alpha=0.5,label=r'Magnetic field computed with $P$ and $\dot{P}$')
plt.legend()
plt.xlabel('Magnetic field (B in Tesla)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_Bfield.pdf')
plt.close()

#Spin-velocity angles histogram
plt.figure(21)
plt.hist(PA,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Number of pulsars')
plt.yscale('log')
plt.savefig('histo_spinvelangle.pdf',dpi=300)
plt.close()

#Spin-velocity angle (old pulsars) histogram
plt.figure(22)
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.xlabel('spin-velocity angle in degrees for pulsars older than 10 Myr')
plt.ylabel('Number of pulsars')
plt.savefig('histo_spinvelangle_old.pdf')
plt.close()

#Spin-velocity angle (young pulsars) histogram
plt.figure(23)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=1,label='Simulation',histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees for pulsars younger than 10 Myr')
plt.ylabel('Number of pulsars')
plt.savefig('histo_spinvelangle_young.pdf')
plt.close()

#Spin-velocity angle (young pulsars and old) histogram
plt.figure(24)
plt.hist(PA_young,bins=20,range=(0,180),edgecolor='blue',color='blue',alpha=1,label='Age < 10Myr',zorder=2,histtype='step')
plt.hist(PA_old,bins=20,range=(0,180),edgecolor='red',color='red',alpha=1,label='Age > 10Myr',histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel('spin-velocity angle in degrees')
plt.ylabel('Number of pulsars')
plt.savefig('histo_spinvelangle_young_and_old.pdf',dpi=300)
plt.close()

#Plot age=f(spin-vel angle)
plt.figure(25)
plt.scatter(PA,log_age,c='red',marker='o',s=5,label='Simulation data')
plt.xlim(0,180)
plt.ylim(0,12)
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Log(age) in yr')
plt.legend()
plt.savefig('spin_vel_age_plot.pdf')
plt.close()

#Initial velocities
plt.figure(26)
plt.hist(v0_norm,bins=20,range=(0,1000),edgecolor='black',color='red',alpha=0.5,label='Simulation')
plt.legend()
plt.xlabel('Birth kick velocities of the pulsars in km/s')
plt.ylabel('Number of pulsars')
plt.savefig('histo_BK_velocities.pdf')
plt.close()

for i in range(len(age)):
    age[i]=age[i]/(365*24*3600) #Get the age in yr

#Plot spin-vel angle=f(nb_orbit)
plt.figure(27)
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
plt.savefig('spinvel_orbit_plot.pdf',dpi=300)
plt.close()

#Plot spin-vel angle=f(v_init)
plt.figure(28)
scatter=plt.scatter(PA,nb_orbit,c=v0,cmap='plasma',marker='o',s=5,label='Simulation data')
plt.ylim(0,10)
plt.xlim(0,180)
plt.xlabel('Spin-velocity angle in degrees')
plt.ylabel('Number of orbits in the Galaxy')
cbar=plt.colorbar(scatter)
cbar.set_label('Initial velocity (km/s)')
plt.legend()
plt.savefig('spinvel_init_vel_plot.pdf',dpi=300)
plt.close()

Fg_flux_log=[]
for i in range(len(Fg_flux)):
    Fg_flux_log.append(np.log10(Fg_flux[i]))

flux_3PC_cano_log=[]
for i in range(len(flux_3PC_cano)):
    flux_3PC_cano_log.append(np.log10(flux_3PC_cano[i]))

#Comparison flux 3PC and simulation (gamma-ray pulsars)
plt.figure(29)
plt.hist(Fg_flux_log,bins=20,range=(-16,-11),edgecolor='red',color='red',alpha=1,label='Simulation',density=True,zorder=3,histtype='step')
plt.hist(flux_3PC_cano_log,bins=20,range=(-16,-11),edgecolor='blue',color='blue',alpha=1,label='3PC data',density=True,zorder=2,histtype='step') #Canonical pop
plt.legend()
plt.yscale('log')
plt.xlabel(r'Gamma-ray flux in log space (W.m$^{-2}$)')
plt.ylabel('Number of pulsars')
plt.savefig('histo_fgamma.pdf',dpi=300)
plt.close()

#Comparison gamma-ray peak separation 3PC and simulation
plt.figure(30)
plt.hist(g_peak_sep_sim,bins=10,range=(0,0.5),edgecolor='red',color='red',alpha=1,label='Simulation',density=True,zorder=3,histtype='step')
plt.hist(g_peak_sep_obs,bins=10,range=(0,0.5),edgecolor='blue',color='blue',alpha=1,label='3PC data',density=True,zorder=2,histtype='step') #Canonical pop
plt.legend()
#plt.yscale('log')
plt.xlabel(r'Gamma-ray peak separation')
plt.ylabel('p.d.f')
plt.savefig('histo_gpeaksep.pdf',dpi=300)
plt.close()

#Rho histogram
plt.figure(31)
plt.hist(rho_rg,bins=20,range=(0,30),edgecolor='blue',color='blue',alpha=0.5,label='Pulsars in both radio and gamma-ray surveys',zorder=3,histtype='step')
plt.hist(rho_r,bins=20,range=(0,30),edgecolor='red',color='red',alpha=0.5,label='Pulsars in radio surveys only',zorder=3,histtype='step')
plt.hist(rho_g,bins=20,range=(0,30),edgecolor='green',color='green',alpha=0.5,label='Pulsars in gamma-ray surveys only',zorder=3,histtype='step')
plt.legend()
plt.xlabel(r'$\rho$ (°)')
plt.ylabel('Number of pulsars')
plt.savefig('rho_sim.pdf',dpi=300)
plt.close()

#Longitude
plt.figure(32)
plt.hist(longitude,bins=20,range=(0,360),edgecolor='red',color='red',alpha=1,label='Simulation',zorder=3,histtype='step')
plt.hist(longitude_ATNF,bins=20,range=(0,360),edgecolor='blue',color='blue',alpha=1,label='ATNF data',zorder=2,histtype='step')
plt.legend()
plt.yscale('log')
plt.xlabel('Longitude in °')
plt.ylabel('Number of pulsars')
plt.savefig('histo_longitude.pdf',dpi=300)
plt.close()

#Zeta=f(chi) all
plt.figure(33)
plt.scatter(alpha_rg,xi_rg,c='blue',marker='o',s=5,label='Pulsars in both radio and gamma-ray surveys')
plt.scatter(alpha_g,xi_g,c='green',marker='x',s=5,label='Pulsars in gamma-ray surveys only')
plt.scatter(alpha_r,xi_r,c='red',marker='s',s=5,label='Pulsars in radio surveys only')
plt.ylabel(r'$\zeta$ in °')
plt.legend(fontsize='small')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.savefig('zeta_chi_all.pdf',dpi=300)
plt.close()

#For w_r with ISM+instrument effects
a_wr,b_wr,r_value,p_value,std_err=linregress(np.log10(P_r_or_rg),np.log10(wr_r_or_rg))
reglinx = np.logspace(np.log10(min(P_r_or_rg)), np.log10(max(P_r_or_rg)), num=len(P_r_or_rg))
regliny = 10**(a_wr * np.log10(reglinx) + b_wr)
residuals_y= np.log10(wr_r_or_rg) - np.log10(regliny)
sigma_y = np.sqrt(np.sum(residuals_y**2) / (len(P_r_or_rg) - 2))
Sxx=np.sum((np.log10(P_r_or_rg) - np.mean(np.log10(P_r_or_rg)))**2)
std_err_b=sigma_y*np.sqrt(1.0/len(P_r_or_rg)+((np.mean(np.log10(P_r_or_rg)))**2/Sxx))
a_wrplus=a_wr+std_err
b_wrplus=b_wr+std_err_b
a_wrminus=a_wr-std_err
b_wrminus=b_wr-std_err_b
reglinyplus = 10**(a_wr * np.log10(reglinx) + b_wrplus)
reglinyminus = 10**(a_wr * np.log10(reglinx) + b_wrminus)

#print(w_geometry_r_or_rg)

#for w_r without ISM+instrument effects
indices=[i for i, val in enumerate(P_r_or_rg)]
P_r_or_rg2=[P_r_or_rg[i] for i in indices]
w_geometry_r_or_rg2=[w_geometry_r_or_rg[i] for i in indices]
a_wr2,b_wr2,r_value2,p_value2,std_err2=linregress(np.log10(P_r_or_rg2),np.log10(w_geometry_r_or_rg2))
reglinx2 = np.logspace(np.log10(min(P_r_or_rg2)), np.log10(max(P_r_or_rg2)), num=len(P_r_or_rg2))
regliny2 = 10**(a_wr2 * np.log10(reglinx2) + b_wr2)
residuals_y2= np.log10(w_geometry_r_or_rg2) - np.log10(regliny2)
sigma_y2 = np.sqrt(np.sum(residuals_y2**2) / (len(P_r_or_rg2) - 2))
Sxx2=np.sum((np.log10(P_r_or_rg2) - np.mean(np.log10(P_r_or_rg2)))**2)
std_err_b2=sigma_y2*np.sqrt(1.0/len(P_r_or_rg2)+((np.mean(np.log10(P_r_or_rg2)))**2/Sxx2))
a_wrplus2=a_wr2+std_err2
b_wrplus2=b_wr2+std_err_b2
a_wrminus2=a_wr2-std_err2
b_wrminus2=b_wr2-std_err_b2
reglinyplus2 = 10**(a_wr2 * np.log10(reglinx2) + b_wrplus2)
reglinyminus2 = 10**(a_wr2 * np.log10(reglinx2) + b_wrminus2)

#W10 Obs
P_obs_combined=P3+P5
w10_combined=w10_r+w10_rg
indices = [i for i, val in enumerate(w10_combined) if val != 0]
w10_combined_filtered = [w10_combined[i] for i in indices]
P_obs_combined_filtered = [P_obs_combined[i] for i in indices]
a_wr3,b_wr3,r_value3,p_value3,std_err3=linregress(np.log10(P_obs_combined_filtered),np.log10(w10_combined_filtered))
reglinx3 = np.logspace(np.log10(min(P_obs_combined_filtered)), np.log10(max(P_obs_combined_filtered)), num=len(P_obs_combined_filtered))
regliny3 = 10**(a_wr3 * np.log10(reglinx3) + b_wr3)
residuals_y3= np.log10(w10_combined_filtered) - np.log10(regliny3)
sigma_y3 = np.sqrt(np.sum(residuals_y3**2) / (len(P_obs_combined_filtered) - 2))
Sxx3=np.sum((np.log10(P_obs_combined_filtered) - np.mean(np.log10(P_obs_combined_filtered)))**2)
std_err_b3=sigma_y3*np.sqrt(1.0/len(P_obs_combined_filtered)+((np.mean(np.log10(P_obs_combined_filtered)))**2/Sxx3))
a_wrplus3=a_wr3+std_err3
b_wrplus3=b_wr3+std_err_b3
a_wrminus3=a_wr3-std_err3
b_wrminus3=b_wr3-std_err_b3
reglinyplus3 = 10**(a_wr3 * np.log10(reglinx3) + b_wrplus3)
reglinyminus3 = 10**(a_wr3 * np.log10(reglinx3) + b_wrminus3)

#print(len(w10_combined_filtered))
#print(len(w_geometry_r_or_rg))

#f(log(wr))=logP
plt.figure(34)
#plt.scatter(P_r_or_rg,wr_r_or_rg,c='green',marker='o',s=5,label='Simulation with ISM and instrumental effect')
plt.scatter(P_r_or_rg2,w_geometry_r_or_rg2,c='red',marker='o',s=5,label='Simulation (geometry only)')
#plt.plot(reglinx,regliny,linestyle='-',label=r'$\log$($w_r$) = (%.2f$\pm%.2f$) $\log(P)$ + (%.2f$\pm$%.2f)' % (a_wr,std_err, b_wr,std_err_b),c='green')
plt.plot(reglinx2,regliny2,linestyle='-',label=r'$\log$($w_r$) = (%.2f$\pm$%.2f) $\log(P)$ + (%.2f$\pm$%.2f) ' % (a_wr2,std_err2,b_wr2,std_err_b2),c='red')
plt.scatter(P3,w10_r,c='blue',marker='o',s=5)
plt.scatter(P5,w10_rg,c='blue',marker='o',s=5,label='ATNF data')
plt.plot(reglinx3,regliny3,linestyle='-',label=r'$\log$($w^{ATNF}_{10}$) = (%.2f$\pm$%.2f) $\log(P)$ + (%.2f$\pm$%.2f)' % (a_wr3,std_err3,b_wr3,std_err_b3),c='blue')
plt.fill_between(reglinx2, reglinyminus2, reglinyplus2, color='red', alpha=0.1)#, label=r'$\pm 1\sigma$ interval')
plt.fill_between(reglinx3, reglinyminus3, reglinyplus3, color='blue', alpha=0.1)
#plt.fill_between(reglinx, reglinyminus, reglinyplus, color='green', alpha=0.1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$w_r$ (°)')
plt.xlabel(r'$P$ (s)')
plt.legend(fontsize='small')
plt.savefig('wr_P_plot.pdf',dpi=300)
plt.close()

#Zeta=f(chi) eliminated
xexample=[i for i in np.arange(0,180,0.1)]
yexample=xexample
plt.figure(35)
#plt.plot(xexample,yexample,label=r'$\rho=\chi+\zeta$')
plt.scatter(alphacheck,xicheck,c='blue',marker='o',s=5)
plt.ylabel(r'$\zeta$ in °')
plt.legend(fontsize='small')
plt.xlabel(r'$\chi$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.savefig('zeta_chi_elim.pdf',dpi=300)
plt.close()

#chi+zeta
alpha_xi2=[]
for i in range(len(alpha_r_or_rg)):
    alpha_xi2.append(alpha_r_or_rg[i]+xi_r_or_rg[i])

plt.figure(36)
plt.scatter(alpha_xi2,rho_r_or_rg,c='blue',marker='o',s=5)
plt.ylabel(r'$\rho$ in °')
plt.legend(fontsize='small')
plt.xlabel(r'$\chi+\zeta$ in °')
plt.xlim(0,90)
plt.ylim(0,90)
plt.axvline(x=20,linestyle='--')
plt.savefig('zeta_chi_elim2.pdf',dpi=300)
plt.close()
