import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits import mplot3d
from scipy.stats import kstest,ks_2samp,linregress
from scipy.stats import mannwhitneyu
from matplotlib.colors import LogNorm
from astropy.table import Table
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib as mpl

def log_format(x, pos):
    return r"$10^{{{}}}$".format(int(x))

def compute_cdf(data, bins):
    hist, edges = np.histogram(data, bins=bins, density=True)
    cdf = np.cumsum(hist) * np.diff(edges)  # CDF normalisée
    return edges[:-1], cdf

#Variable initialization
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
Pa,P_dota,da,za,xa,ya,agea,E_dota=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue, before the unknown values are ruled out 
Pa_gamma,P_dota_gamma,da_gamma,za_gamma,xa_gamma,ya_gamma,agea_gamma,E_dota_gamma=[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue of the gamma pulsars only
P2,P_dot2,d2,z2,x2,y2,age2,E_dot2,latitude2,longitude2=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used
P3,P_dot3,d3,z3,x3,y3,age3,E_dot3,latitude3,longitude3=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio)
P4,P_dot4,d4,z4,x4,y4,age4,E_dot4,latitude4,longitude4=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (gamma)
P5,P_dot5,d5,z5,x5,y5,age5,E_dot5,latitude5,longitude5=[],[],[],[],[],[],[],[],[],[] #Refers to the data of the ATNF catalogue that is used (radio-gamma)
P_x,P_dot_x,x_x,y_x,age_x,error_x,type_pulsar_x,distance_x,latitude_x,longitude_x,cos_alpha0_x,cos_alpha_x,Bf_x,z_x,vx_x,vy_x,vz_x,vx0_x,vy0_x,vz0_x,PA_x=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data (x-ray)
log_age,log_age2,log_P,log_Pdot,log_P2,log_Pdot2=[],[],[],[],[],[] #Refers to the quantities we need in log scale
RAD = 180/np.pi
test_l=[]
B4,B5,B3=[],[],[]
Ba,B2,Ba_gamma=[],[],[]
B0=3.2e15 #constant to compute the surface magnetic field of the observations
P_selected,Pdot_selected=[],[]
#B_init=[]

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

#Simulation data
P,P_dot,x,y,age,error,type_pulsar,distance,latitude,longitude,cos_alpha0,cos_alpha,Bf,z,vx,vy,vz,vx0,vy0,vz0,PA=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[] #Refers to the simulation data
var,var2='',''
reg_1=re.compile("-*.{12}[|]{1}")
reg_2=re.compile("[|]{1}.{1}[|]{1}")
reg_3=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")

with open("x_file.txt","r") as f:
    data_x=re.findall(reg_1,f.read())

with open("x_file.txt","r") as f:
    data_type_x=re.findall(reg_2,f.read())

with open("x_file2.txt","r") as f:
    data2_x=re.findall(reg_3,f.read())

with open("PF_check.txt","r") as f:
    data_PF=re.findall(reg_3,f.read())

Fxmax,Fxmin,PF=[],[],[]
for i in range(int(len(data_PF)/3)):
    PF.append(float(data_PF[3*i+0]))
    Fxmax.append(float(data_PF[3*i+1]))
    Fxmin.append(float(data_PF[3*i+2]))

for i in range(len(data_type_x)):
    var2+=data_type_x[i][1]
    type_pulsar_x+=[int(var2)]
    var2=''

for i in range(len(data_x)):
    for j in range(len(data_x[i])-1):
        var+=data_x[i][j]
    data_x[i]=float(var)
    var=''

for i in range(int(len(data_x)/20)):
    P_x+=[data_x[20*i]]    #Period of rotation of the pulsar in seconds
    P_dot_x+=[data_x[20*i+1]]  #Period derivative of rotation of the pulsar, no units
    x_x+=[data_x[20*i+2]]#Position on the x-absciss relative to the sun in the galactocentric frame in kpc
    y_x+=[data_x[20*i+3]]#Position on the y-absciss relative to the sun in the galactocentric frame in kpc
    age_x+=[data_x[20*i+4]] #Age of the pulsar in seconds
    error_x+=[data_x[20*i+5]] #error on the positions of the pulsars (error computed with the energy)
    distance_x+=[data_x[20*i+6]] #Distance to the galactic center in kpc
    longitude_x+=[data_x[20*i+7]] #galactic latitude in degrees
    latitude_x+=[data_x[20*i+8]] #galactic longitude in degrees
    cos_alpha0_x+=[data_x[20*i+9]] #cosinus of the initial inclination angle
    cos_alpha_x+=[np.cos(data_x[20*i+10])] #cosinus of the inclination angle after all the evolution
    Bf_x+=[data_x[20*i+11]] #Magnetic field of the pulsar today, in tesla
    z_x+=[data_x[20*i+12]] #Position on the z-absciss relative to the sun in the galactocentric frame in kpc
    vx_x+=[data_x[20*i+13]] #Velocity on the x-absciss of the pulsar today, in km/s
    vy_x+=[data_x[20*i+14]] #Velocity on the y-absciss of the pulsar today, in km/s
    vz_x+=[data_x[20*i+15]] #Velocity on the z-absciss of the pulsar today, in km/s
    vx0_x+=[data_x[20*i+16]] #Velocity on the x-absciss of the pulsar initially, in km/s
    vy0_x+=[data_x[20*i+17]] #Velocity on the y-absciss of the pulsar initially, in km/s
    vz0_x+=[data_x[20*i+18]] #Velocity on the z-absciss of the pulsar initially, in km/s
    PA_x+=[(180/np.pi)*np.arccos(data_x[20*i+19])]

T_n,T_s,r_hn,cos_in,F_x,xi,r_hs=[],[],[],[],[],[],[]
cos_is,theta_n,phi_n,theta_s,phi_s,mu_hot_ang1,mu_hot_ang2=[],[],[],[],[],[],[]
for i in range(int(len(data2_x)/14)):
    cos_in.append(float(data2_x[14*i])) #cos of the angle between magnetic axis and the line of sight (north)
    T_n.append(float(data2_x[14*i+1])) #Temperature of the hot spot in K
    r_hn.append(float(data2_x[14*i+2])) #Radius of the hot spot in meter
    F_x.append(float(data2_x[14*i+3])) #Flux in thermal X-ray in W.m^-2
    xi.append(float(data2_x[14*i+4])*180/np.pi) #Viewing angle in degree
    cos_is.append(float(data2_x[14*i+5])) #cos of the angle between magnetic axis and the line of sight (south)
    theta_n.append(float(data2_x[14*i+6])) #Position in spherical coordinate in rad (theta and north)
    phi_n.append(float(data2_x[14*i+7])) #Position in spherical coordinate in rad (phi and north)
    theta_s.append(float(data2_x[14*i+8])) #Position in spherical coordinate in rad (theta and south)
    phi_s.append(float(data2_x[14*i+9])) #Position in spherical coordinate in rad (phi and south)
    mu_hot_ang1.append(float(data2_x[14*i+10])*180/np.pi) #Angle in degree between magnetic axis and hotspot (north)
    mu_hot_ang2.append(float(data2_x[14*i+11])*180/np.pi) #Angle in degree between magnetic axis and hotspot (south)
    r_hs.append(float(data2_x[14*i+12])) #Radius of the hot spot in meter (south)
    T_s.append(float(data2_x[14*i+13])) #Temperature of the hot spot in K (south)

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
    Lx_abs_redshift.append(F_x[i]*4*np.pi*(distance_x[i]*3.086e19)**2)

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

with open("init_P_B.txt","r") as f:
    init_P_B=re.findall(reg_5,f.read())

with open("Lgamma.txt","r") as f:
    Lgamma_saved=re.findall(reg_4,f.read())

with open("PPdot_alpha_all.txt","r") as f:
    all_ppdot_alpha=re.findall(reg_5,f.read())

df=pd.read_excel('3PC_Catalog_20230803.xls')
data_3PC=Table.from_pandas(df)

#Data Alex
#alpha_tab=pd.read_csv("alpha_vs_aligntime_alex.csv")
#rho_tab=pd.read_csv("rho_vs_P_alex.csv")

#tc_alex_sample=((alpha_tab['Period(s)'])/(2*alpha_tab['Pdot']))/(365*24*3600)

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

#Get the Lgamma data
Lgamma_saved2=[]
for i in range(int(len(Lgamma_saved))):
    Lgamma_saved2.append(float(Lgamma_saved[i])*1e7) #In erg/s

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

#Get the data from the file all_ppdot_alpha.txt
P_sim_all,P_dot_sim_all,cos_alpha_all,cos_alpha0_all=[],[],[],[]
for i in range(int(len(all_ppdot_alpha)/4)):
    P_sim_all.append(float(all_ppdot_alpha[4*i+0]))
    P_dot_sim_all.append(float(all_ppdot_alpha[4*i+1]))
    cos_alpha_all.append(float(np.cos(float(all_ppdot_alpha[4*i+2]))))
    cos_alpha0_all.append(float(all_ppdot_alpha[4*i+3]))

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

#Get initial P and B from sim
P0_sim,B0_sim=[],[]
for i in range(int(len(init_P_B)/2)):
    B0_sim.append(float(init_P_B[2*i]))
    P0_sim.append(float(init_P_B[2*i+1]))

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

with open("info_supp_obs.txt","w") as f:
    f.write(f'{len(P3)}\n')
print(f'Number of observed radio pulsars (FAST GPPS + PMPS): {len(P3)}')
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
R_NS=10000
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

#List of magnetic obliquity angle
alpha_all,tau_MHD_align0,tau_MHD_align,alpha_all0,t0_tau,tc_tau=[],[],[],[],[],[]
for i in range(len(cos_alpha)):
    al=min(180*np.arccos(cos_alpha[i])/np.pi,180-180*np.arccos(cos_alpha[i])/np.pi)
    al0=min(180*np.arccos(cos_alpha0[i])/np.pi,180-180*np.arccos(cos_alpha0[i])/np.pi)
    alpha_all.append(al)
    alpha_all0.append(al0)
    tau_MHD_align0.append(np.log10(((Inertia*mu_0*c_light**3*P0_sim[i]**2*np.sin(al0*np.pi/180)**2)/(16*np.pi**3*R_NS**6*B0_sim[i]**2*np.cos(al0*np.pi/180)**4))/(365*24*3600))) #in yr + logscale
    tau_MHD_align.append(np.log10(((Inertia*mu_0*c_light**3*P[i]**2*np.sin(al*np.pi/180)**2)/(16*np.pi**3*R_NS**6*Bf[i]**2*np.cos(al*np.pi/180)**4))/(365*24*3600))) #in yr + logscale

#tau_MHD_alex_sample=((Inertia*mu_0*c_light**3*alpha_tab['Period(s)']**2*np.sin(alpha_tab['Alpha(deg)']*np.pi/180)**2)/(16*np.pi**3*R_NS**6*(((Inertia*mu_0*c_light**3)/(16*np.pi**3*R_NS**6*(1+np.sin(alpha_tab['Alpha(deg)']*np.pi/180)**2)))*alpha_tab['Period(s)']*alpha_tab['Pdot'])**1*np.cos(alpha_tab['Alpha(deg)']*np.pi/180)**4))/(365*24*3600)
#ratio_tc_tau_MHD_alex=tc_alex_sample/tau_MHD_alex_sample
#ratio_tau_MHD_tc_alex=tau_MHD_alex_sample/tc_alex_sample

for i in range(len(age)):
    tc_tau.append((charac_age[i])/(10**tau_MHD_align[i])) #Ratio age/tau_MHD_align

tau_tc=[]
for i in range(len(age)):
    tau_tc.append((10**tau_MHD_align[i])/(charac_age[i])) #Ratio tau_MHD_align/age_charac

#Simulated sample
countless1e5,count1e5,count1e6,count1e7,count1e8,count1e9,count1e10=0,0,0,0,0,0,0
for i in range(len(tau_MHD_align)):
    if tau_MHD_align[i]<5:
        countless1e5+=1
    elif tau_MHD_align[i]>=5 and tau_MHD_align[i]<6:
        count1e5+=1
    elif tau_MHD_align[i]>=6 and tau_MHD_align[i]<7:
        count1e6+=1
    elif tau_MHD_align[i]>=7 and tau_MHD_align[i]<8:
        count1e7+=1
    elif tau_MHD_align[i]>=8 and tau_MHD_align[i]<9:
        count1e8+=1
    elif tau_MHD_align[i]>=9 and tau_MHD_align[i]<10:
        count1e9+=1
    elif tau_MHD_align[i]>=10:
        count1e10+=1

print(f'----------SIMULATION-----------')
print(f'Proportion of pulsars with an alignment timescale below 1e5 yr: {countless1e5/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale between 1e5 and 1e6 yr: {count1e5/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale between 1e6 and 1e7 yr: {count1e6/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale between 1e7 and 1e8 yr: {count1e7/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale between 1e8 and 1e9 yr: {count1e8/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale between 1e9 and 1e10 yr: {count1e9/len(tau_MHD_align)}')
print(f'Proportion of pulsars with an alignment timescale above 1e10 yr: {count1e10/len(tau_MHD_align)}')

#Observed sample
#countless1e5,count1e5,count1e6,count1e7,count1e8,count1e9,count1e10=0,0,0,0,0,0,0
#tau_MHD_alex_sample=np.array(tau_MHD_alex_sample)
#for i in range(len(tau_MHD_alex_sample)):
#    if tau_MHD_alex_sample[i]<1e5:
#        countless1e5+=1
#    elif tau_MHD_alex_sample[i]>=1e5 and tau_MHD_alex_sample[i]<1e6:
#        count1e5+=1
#    elif tau_MHD_alex_sample[i]>=1e6 and tau_MHD_alex_sample[i]<1e7:
#        count1e6+=1
#    elif tau_MHD_alex_sample[i]>=1e7 and tau_MHD_alex_sample[i]<1e8:
#        count1e7+=1
#    elif tau_MHD_alex_sample[i]>=1e8 and tau_MHD_alex_sample[i]<1e9:
#        count1e8+=1
#    elif tau_MHD_alex_sample[i]>=1e9 and tau_MHD_alex_sample[i]<1e10:
#        count1e9+=1
#    elif tau_MHD_alex_sample[i]>=1e10:
#        count1e10+=1

#print(f'----------OBSERVATION-----------')
#print(f'Proportion of pulsars with an alignment timescale below 1e5 yr: {countless1e5/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale between 1e5 and 1e6 yr: {count1e5/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale between 1e6 and 1e7 yr: {count1e6/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale between 1e7 and 1e8 yr: {count1e7/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale between 1e8 and 1e9 yr: {count1e8/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale between 1e9 and 1e10 yr: {count1e9/len(tau_MHD_alex_sample)}')
#print(f'Proportion of pulsars with an alignment timescale above 1e10 yr: {count1e10/len(tau_MHD_alex_sample)}')

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
#KS_test=kstest(P_dot,P_dota)
#test_stat_all1=KS_test.statistic
#p_value_all1=KS_test.pvalue
#print("----ALL THE DATA----\n")
#print(f"d_value of Pdot KS test = {test_stat_all1}")
#print(f"p_value of Pdot={p_value_all1}")

#KS_test=kstest(P,Pa)
#test_stat_all2=KS_test.statistic
#p_value_all2=KS_test.pvalue
#print(f"d_value of P KS test = {test_stat_all2}")
#print(f"p_value of P={p_value_all2}")

#KS test 1D python (gamma-ray only population)
#KS_test=kstest(P_dot_gamma,data_3PC_filtered['P1'])
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print("----GAMMA ONLY PULSARS----\n")
#print(f"d_value of Pdot KS test for the gamma only pulsars= {test_stat}")
#print(f"p_value of Pdot for the gamma only pulsars={p_value}")

#KS_test=kstest(P_gamma,data_3PC_filtered['P0'])
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print(f"d_value of P KS test for the gamma only pulsars= {test_stat}")
#print(f"p_value of P for the gamma only pulsars={p_value}")

#KS test 1D python (radio/gamma-ray population)
#KS_test=kstest(P_dot_radio_gamma,P_dot5)
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print("----ALL RADIO/GAMMA PULSARS----\n")
#print(f"d_value of Pdot KS test for the radio/gamma pulsars= {test_stat}")
#print(f"p_value of Pdot for the radio/gamma pulsars={p_value}")

#KS_test=kstest(P_radio_gamma,P5)
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print(f"d_value of P KS test for the radio/gamma pulsars= {test_stat}")
#print(f"p_value of P for the radio/gamma pulsars={p_value}")

#KS test 1D python (radio population)
#KS_test=kstest(P_dot_radio,P_dot3)
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print("----ALL RADIO ONLY PULSARS----\n")
#print(f"d_value of Pdot KS test for the radio only pulsars= {test_stat}")
#print(f"p_value of Pdot for the radio only pulsars={p_value}")

#KS_test=kstest(P_radio,P3)
#test_stat=KS_test.statistic
#p_value=KS_test.pvalue
#print(f"d_value of P KS test for the radio only pulsars= {test_stat}")
#print(f"p_value of P for the radio only pulsars={p_value}")

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

#P_obs_all=np.concatenate(np.array(Pa,dtype=float),np.array(data_X2['P'],dtype=float))
P_obs_all = np.concatenate([
    np.array(Pa, dtype=float),
    np.array(data_X2['P'], dtype=float)
])
#P_dot_obs_all=np.concatenate(np.array(P_dota,dtype=float),np.array(data_X2['Pdot'],dtype=float))
P_dot_obs_all = np.concatenate([
    np.array(P_dota, dtype=float),
    np.array(data_X2['Pdot'], dtype=float)
])

#Computation of the CDF of P and Pdot, obs and sim all
bins_P = np.linspace(-2, np.log10(3e1), 33)
bins_Pdot = np.linspace(-19, -10, 33)

sim_P_cdf_x, sim_P_cdf_y = compute_cdf(np.log10(P_sim_all), bins_P)
obs_P_cdf_x, obs_P_cdf_y = compute_cdf(np.log10(P_obs_all), bins_P)

sim_Pdot_cdf_x, sim_Pdot_cdf_y = compute_cdf(np.log10(P_dot_sim_all), bins_Pdot)
obs_Pdot_cdf_x, obs_Pdot_cdf_y = compute_cdf(np.log10(P_dot_obs_all), bins_Pdot)

sim_P_cdf_x = 10**sim_P_cdf_x
obs_P_cdf_x = 10**obs_P_cdf_x

sim_Pdot_cdf_x = 10**sim_Pdot_cdf_x
obs_Pdot_cdf_x = 10**obs_Pdot_cdf_x

bins_P = np.logspace(-2, np.log10(3e1), 33)
bins_Pdot = np.logspace(-19, -10, 33)

#Plots
#2D histogram P-Pdot all
plt.figure(0)
weights = np.ones_like(P_sim_all) / len(P_sim_all)
hist=plt.hist2d(P_sim_all,P_dot_sim_all,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
H,X,Y,_=plt.hist2d(P_sim_all,P_dot_sim_all,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Fraction of occurrences from simulations')
plt.xlabel(r'Spin Period [s]')
plt.ylabel(r'Period Derivative [s/s]')
plt.gca().xaxis.set_major_formatter(FuncFormatter(log_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(log_format))
plt.savefig(f'inference_plots/2D_P_Pdot_plot_sim.pdf',dpi=300)
plt.close()

#Contour map P-Pdot simulation all
plt.figure(1,figsize=(8,8))
ax_main = plt.axes([0.1, 0.1, 0.65, 0.65])
X_grid, Y_grid = np.meshgrid(X[:-1], Y[:-1])
contour=ax_main.contour(X_grid, Y_grid, H.T, levels=10, cmap='viridis',zorder=2)  # Transposer H car les axes sont inversés
#ax_main.xaxis.set_major_formatter(FuncFormatter(log_format))
#ax_main.set_xticks([-2, -1, 0, 1, 2])
#ax_main.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',r'$10^{1}$',r'$10^{2}$'])
#ax_main.yaxis.set_major_formatter(FuncFormatter(log_format))
#cbar=plt.colorbar(contour)
#cbar.set_label('Fraction of occurrences from simulations')
ax_main.scatter(P_obs_all,P_dot_obs_all,color='red',marker="*",s=15,zorder=1)
ax_main.set_xlabel(r'Spin Period [s]')
ax_main.set_ylabel(r'Period Derivative [s/s]')
ax_main.set_xlim(1e-2,3e1)
ax_main.set_ylim(1e-19,1e-10)
ax_main.set_yscale('log')
ax_main.set_xscale('log')
#ax_main.legend()
#ax_main.grid(alpha=0.5, linestyle='-')
#Bottom CDF (for P)
ax_top = plt.axes([0.1, 0.78, 0.65, 0.1], sharex=ax_main)  # Aligné avec le contour
ax_top.plot(sim_P_cdf_x, sim_P_cdf_y, label='Simulated', color='green')
ax_top.plot(obs_P_cdf_x, obs_P_cdf_y, label='Observed', color='red')
ax_top.legend()
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.set_ylabel("CDF")
#ax_top.set_yscale('log')
#ax_top.set_xscale('log')
ax_top.grid(alpha=0.5, linestyle='-')
#Right CDF (for Pdot)
ax_right = plt.axes([0.77, 0.1, 0.1, 0.65], sharey=ax_main)  # Aligné avec le contour
ax_right.plot(sim_Pdot_cdf_y, sim_Pdot_cdf_x, color='green')  # Transposée
ax_right.plot(obs_Pdot_cdf_y, obs_Pdot_cdf_x, color='red')  # Transposée
#ax_right.set_yscale('log')
#ax_right.set_xscale('log')
#ax_right.legend()
ax_right.tick_params(axis="y", labelleft=False)
ax_right.set_xlabel("CDF")
ax_right.grid(alpha=0.5, linestyle='-')
#Plot the Edot lines
ax_main.plot(P_line,P_dot_Edot1e22W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[1800]*0.85,P_dot_Edot1e22W[1800]*1.05,r'$\dot{E} = 10^{22}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[800]*0.85,P_dot_Edot1e25W[800]*1.05,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[77]*0.85,P_dot_Edot1e28W[77]*1.05,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[7]*0.85,P_dot_Edot1e31W[7]*1.05,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
ax_main.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[80]*0.6,P_dot_B1e6[80]*0.95,r'$B =10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[1500]*0.6,P_dot_B5e6[1500]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2000]*0.6,P_dot_B1e7[2000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2300]*0.6,P_dot_B1e8[2300]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e9,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2500]*0.6,P_dot_B1e9[2500]*0.95,r'$B= 10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B4_4e9,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[2500]*0.6,P_dot_B4_4e9[2500]*0.95,r'$B= 4.4\times10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.savefig(f'inference_plots/contour_map_ppdot_all.pdf',dpi=300)

#Computation of the CDF of P and Pdot, obs and sim gamma
bins_P = np.linspace(-2, np.log10(3e1), 33)
bins_Pdot = np.linspace(-19, -10, 33)

sim_P_cdf_x, sim_P_cdf_y = compute_cdf(np.log10(P_gamma), bins_P)
obs_P_cdf_x, obs_P_cdf_y = compute_cdf(np.log10(data_3PC_filtered['P0']), bins_P)

sim_Pdot_cdf_x, sim_Pdot_cdf_y = compute_cdf(np.log10(P_dot_gamma), bins_Pdot)
obs_Pdot_cdf_x, obs_Pdot_cdf_y = compute_cdf(np.log10(data_3PC_filtered['P1']), bins_Pdot)

sim_P_cdf_x = 10**sim_P_cdf_x
obs_P_cdf_x = 10**obs_P_cdf_x

sim_Pdot_cdf_x = 10**sim_Pdot_cdf_x
obs_Pdot_cdf_x = 10**obs_Pdot_cdf_x

bins_P = np.logspace(-2, np.log10(3e1), 33)
bins_Pdot = np.logspace(-19, -10, 33)

#2D histogram P-Pdot gamma
plt.figure(2)
weights = np.ones_like(P_gamma) / len(P_gamma)
hist=plt.hist2d(P_gamma,P_dot_gamma,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
H,X,Y,_=plt.hist2d(P_gamma,P_dot_gamma,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Fraction of occurrences from simulations')
plt.xlabel(r'Spin Period [s]')
plt.ylabel(r'Period Derivative [s/s]')
plt.gca().xaxis.set_major_formatter(FuncFormatter(log_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(log_format))
plt.savefig(f'inference_plots/2D_P_Pdot_plot_gamma_sim.pdf',dpi=300)
plt.close()

#Contour map P-Pdot simulation gamma
plt.figure(3,figsize=(8,8))
ax_main = plt.axes([0.1, 0.1, 0.65, 0.65])
X_grid, Y_grid = np.meshgrid(X[:-1], Y[:-1])
contour=ax_main.contour(X_grid, Y_grid, H.T, levels=10, cmap='viridis')  # Transposer H car les axes sont inversés
#ax_main.xaxis.set_major_formatter(FuncFormatter(log_format))
#ax_main.set_xticks([-2, -1, 0, 1, 2])
#ax_main.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',r'$10^{1}$',r'$10^{2}$'])
#ax_main.yaxis.set_major_formatter(FuncFormatter(log_format))
#cbar=plt.colorbar(contour)
#cbar.set_label('Fraction of occurrences from simulations')
ax_main.scatter(data_3PC_filtered['P0'],data_3PC_filtered['P1'],color='red',marker="*",s=15,zorder=2)
ax_main.set_xlabel(r'Spin Period [s]')
ax_main.set_ylabel(r'Period Derivative [s/s]')
ax_main.set_xlim(1e-2,3e1)
ax_main.set_ylim(1e-19,1e-10)
ax_main.set_yscale('log')
ax_main.set_xscale('log')
#ax_main.legend()
#ax_main.grid(alpha=0.5, linestyle='-')
#Bottom CDF (for P)
ax_top = plt.axes([0.1, 0.78, 0.65, 0.1], sharex=ax_main)  # Aligné avec le contour
ax_top.plot(sim_P_cdf_x, sim_P_cdf_y, label='Simulated', color='green')
ax_top.plot(obs_P_cdf_x, obs_P_cdf_y, label='Observed', color='red')
ax_top.legend()
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.set_ylabel("CDF")
ax_top.grid(alpha=0.5, linestyle='-')
#Right CDF (for Pdot)
ax_right = plt.axes([0.77, 0.1, 0.1, 0.65], sharey=ax_main)  # Aligné avec le contour
ax_right.plot(sim_Pdot_cdf_y, sim_Pdot_cdf_x, color='green')  # Transposée
ax_right.plot(obs_Pdot_cdf_y, obs_Pdot_cdf_x, color='red')  # Transposée
#ax_right.legend()
ax_right.tick_params(axis="y", labelleft=False)
ax_right.set_xlabel("CDF")
ax_right.grid(alpha=0.5, linestyle='-')
#Plot the Edot lines
ax_main.plot(P_line,P_dot_Edot1e22W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[1800]*0.85,P_dot_Edot1e22W[1800]*1.05,r'$\dot{E} = 10^{22}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[800]*0.85,P_dot_Edot1e25W[800]*1.05,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[77]*0.85,P_dot_Edot1e28W[77]*1.05,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[7]*0.85,P_dot_Edot1e31W[7]*1.05,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
ax_main.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[80]*0.6,P_dot_B1e6[80]*0.95,r'$B =10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[1500]*0.6,P_dot_B5e6[1500]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2000]*0.6,P_dot_B1e7[2000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2300]*0.6,P_dot_B1e8[2300]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e9,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2500]*0.6,P_dot_B1e9[2500]*0.95,r'$B= 10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B4_4e9,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[2500]*0.6,P_dot_B4_4e9[2500]*0.95,r'$B= 4.4\times10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.savefig(f'inference_plots/contour_map_ppdot_gamma.pdf',dpi=300)

#Computation of the CDF of P and Pdot, obs and sim radio
bins_P = np.linspace(-2, np.log10(3e1), 33)
bins_Pdot = np.linspace(-19, -10, 33)

sim_P_cdf_x, sim_P_cdf_y = compute_cdf(np.log10(P_radio), bins_P)
obs_P_cdf_x, obs_P_cdf_y = compute_cdf(np.log10(P3), bins_P)

sim_Pdot_cdf_x, sim_Pdot_cdf_y = compute_cdf(np.log10(P_dot_radio), bins_Pdot)
obs_Pdot_cdf_x, obs_Pdot_cdf_y = compute_cdf(np.log10(P_dot3), bins_Pdot)

sim_P_cdf_x = 10**sim_P_cdf_x
obs_P_cdf_x = 10**obs_P_cdf_x

sim_Pdot_cdf_x = 10**sim_Pdot_cdf_x
obs_Pdot_cdf_x = 10**obs_Pdot_cdf_x

bins_P = np.logspace(-2, np.log10(3e1), 33)
bins_Pdot = np.logspace(-19, -10, 33)

#2D histogram P-Pdot radio
plt.figure(4)
weights = np.ones_like(P_radio) / len(P_radio)
hist=plt.hist2d(P_radio,P_dot_radio,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
H,X,Y,_=plt.hist2d(P_radio,P_dot_radio,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Fraction of occurrences from simulations')
plt.xlabel(r'Spin Period [s]')
plt.ylabel(r'Period Derivative [s/s]')
plt.gca().xaxis.set_major_formatter(FuncFormatter(log_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(log_format))
plt.savefig(f'inference_plots/2D_P_Pdot_plot_radio_sim.pdf',dpi=300)
plt.close()

#Contour map P-Pdot simulation gamma
plt.figure(5,figsize=(8,8))
ax_main = plt.axes([0.1, 0.1, 0.65, 0.65])
X_grid, Y_grid = np.meshgrid(X[:-1], Y[:-1])
contour=ax_main.contour(X_grid, Y_grid, H.T, levels=10, cmap='viridis')  # Transposer H car les axes sont inversés
#ax_main.xaxis.set_major_formatter(FuncFormatter(log_format))
#ax_main.set_xticks([-2, -1, 0, 1, 2])
#ax_main.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',r'$10^{1}$',r'$10^{2}$'])
#ax_main.yaxis.set_major_formatter(FuncFormatter(log_format))
#cbar=plt.colorbar(contour)
#cbar.set_label('Fraction of occurrences from simulations')
ax_main.scatter(P3,P_dot3,color='red',marker="*",s=15,zorder=2)
ax_main.set_xlabel(r'Spin Period [s]')
ax_main.set_ylabel(r'Period Derivative [s/s]')
ax_main.set_xlim(1e-2,3e1)
ax_main.set_ylim(1e-19,1e-10)
ax_main.set_yscale('log')
ax_main.set_xscale('log')
#ax_main.legend()
#ax_main.grid(alpha=0.5, linestyle='-')
#Bottom CDF (for P)
ax_top = plt.axes([0.1, 0.78, 0.65, 0.1], sharex=ax_main)  # Aligné avec le contour
ax_top.plot(sim_P_cdf_x, sim_P_cdf_y, label='Simulated', color='green')
ax_top.plot(obs_P_cdf_x, obs_P_cdf_y, label='Observed', color='red')
ax_top.legend()
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.set_ylabel("CDF")
ax_top.grid(alpha=0.5, linestyle='-')
#Right CDF (for Pdot)
ax_right = plt.axes([0.77, 0.1, 0.1, 0.65], sharey=ax_main)  # Aligné avec le contour
ax_right.plot(sim_Pdot_cdf_y, sim_Pdot_cdf_x, color='green')  # Transposée
ax_right.plot(obs_Pdot_cdf_y, obs_Pdot_cdf_x, color='red')  # Transposée
#ax_right.legend()
ax_right.tick_params(axis="y", labelleft=False)
ax_right.set_xlabel("CDF")
ax_right.grid(alpha=0.5, linestyle='-')
#Plot the Edot lines
ax_main.plot(P_line,P_dot_Edot1e22W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[1800]*0.85,P_dot_Edot1e22W[1800]*1.05,r'$\dot{E} = 10^{22}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[800]*0.85,P_dot_Edot1e25W[800]*1.05,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[77]*0.85,P_dot_Edot1e28W[77]*1.05,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[7]*0.85,P_dot_Edot1e31W[7]*1.05,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
ax_main.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[80]*0.6,P_dot_B1e6[80]*0.95,r'$B =10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[1500]*0.6,P_dot_B5e6[1500]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2000]*0.6,P_dot_B1e7[2000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2300]*0.6,P_dot_B1e8[2300]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e9,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2500]*0.6,P_dot_B1e9[2500]*0.95,r'$B= 10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B4_4e9,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[2500]*0.6,P_dot_B4_4e9[2500]*0.95,r'$B= 4.4\times10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.savefig(f'inference_plots/contour_map_ppdot_radio.pdf',dpi=300)

#Computation of the CDF of P and Pdot, obs and sim thermal x-rays
bins_P = np.linspace(-2, np.log10(3e1), 33)
bins_Pdot = np.linspace(-19, -10, 33)

sim_P_cdf_x, sim_P_cdf_y = compute_cdf(np.log10(P_x), bins_P)
obs_P_cdf_x, obs_P_cdf_y = compute_cdf(np.log10(data_X2['P']), bins_P)

sim_Pdot_cdf_x, sim_Pdot_cdf_y = compute_cdf(np.log10(P_dot_x), bins_Pdot)
obs_Pdot_cdf_x, obs_Pdot_cdf_y = compute_cdf(np.log10(data_X2['Pdot']), bins_Pdot)

sim_P_cdf_x = 10**sim_P_cdf_x
obs_P_cdf_x = 10**obs_P_cdf_x

sim_Pdot_cdf_x = 10**sim_Pdot_cdf_x
obs_Pdot_cdf_x = 10**obs_Pdot_cdf_x

bins_P = np.logspace(-2, np.log10(3e1), 33)
bins_Pdot = np.logspace(-19, -10, 33)

#2D histogram P-Pdot x-ray
plt.figure(6)
weights = np.ones_like(P_x) / len(P_x)
hist=plt.hist2d(P_x,P_dot_x,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
H,X,Y,_=plt.hist2d(P_x,P_dot_x,bins=(bins_P,bins_Pdot),cmap='viridis',weights=weights)
cbar=plt.colorbar(hist[3])
cbar.set_label('Fraction of occurrences from simulations')
plt.xlabel(r'Spin Period [s]')
plt.ylabel(r'Period Derivative [s/s]')
plt.gca().xaxis.set_major_formatter(FuncFormatter(log_format))
plt.gca().yaxis.set_major_formatter(FuncFormatter(log_format))
plt.savefig(f'inference_plots/2D_P_Pdot_plot_x_sim.pdf',dpi=300)
plt.close()

#Contour map P-Pdot simulation thermal x-ray
plt.figure(7,figsize=(8,8))
ax_main = plt.axes([0.1, 0.1, 0.65, 0.65])
X_grid, Y_grid = np.meshgrid(X[:-1], Y[:-1])
contour=ax_main.contour(X_grid, Y_grid, H.T, levels=10, cmap='viridis')  # Transposer H car les axes sont inversés
#ax_main.xaxis.set_major_formatter(FuncFormatter(log_format))
#ax_main.set_xticks([-2, -1, 0, 1, 2])
#ax_main.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$',r'$10^{1}$',r'$10^{2}$'])
#ax_main.yaxis.set_major_formatter(FuncFormatter(log_format))
#cbar=plt.colorbar(contour)
#cbar.set_label('Fraction of occurrences from simulations')
ax_main.scatter(data_X2['P'],data_X2['Pdot'],color='red',marker="*",s=15,zorder=2)
ax_main.set_xlabel(r'Spin Period [s]')
ax_main.set_ylabel(r'Period Derivative [s/s]')
ax_main.set_xlim(1e-2,3e1)
ax_main.set_ylim(1e-19,1e-10)
ax_main.set_yscale('log')
ax_main.set_xscale('log')
#ax_main.legend()
#ax_main.grid(alpha=0.5, linestyle='-')
#Bottom CDF (for P)
ax_top = plt.axes([0.1, 0.78, 0.65, 0.1], sharex=ax_main)  # Aligné avec le contour
ax_top.plot(sim_P_cdf_x, sim_P_cdf_y, label='Simulated', color='green')
ax_top.plot(obs_P_cdf_x, obs_P_cdf_y, label='Observed', color='red')
ax_top.legend()
ax_top.tick_params(axis="x", labelbottom=False)
ax_top.set_ylabel("CDF")
ax_top.grid(alpha=0.5, linestyle='-')
#Right CDF (for Pdot)
ax_right = plt.axes([0.77, 0.1, 0.1, 0.65], sharey=ax_main)  # Aligné avec le contour
ax_right.plot(sim_Pdot_cdf_y, sim_Pdot_cdf_x, color='green')  # Transposée
ax_right.plot(obs_Pdot_cdf_y, obs_Pdot_cdf_x, color='red')  # Transposée
#ax_right.legend()
ax_right.tick_params(axis="y", labelleft=False)
ax_right.set_xlabel("CDF")
ax_right.grid(alpha=0.5, linestyle='-')
#Plot the Edot lines
ax_main.plot(P_line,P_dot_Edot1e22W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[1800]*0.85,P_dot_Edot1e22W[1800]*1.05,r'$\dot{E} = 10^{22}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e25W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[800]*0.85,P_dot_Edot1e25W[800]*1.05,r'$\dot{E} = 10^{25}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e28W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[77]*0.85,P_dot_Edot1e28W[77]*1.05,r'$\dot{E} = 10^{28}$W',fontsize=7,color='orange',rotation=45,zorder=3)
ax_main.plot(P_line,P_dot_Edot1e31W,linestyle='dotted',c='orange',zorder=0)
ax_main.text(P_line[7]*0.85,P_dot_Edot1e31W[7]*1.05,r'$\dot{E} = 10^{31}$W',fontsize=7,color='orange',rotation=45,zorder=3)
#Plot the B lines
ax_main.plot(P_line,P_dot_B1e6,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[80]*0.6,P_dot_B1e6[80]*0.95,r'$B =10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B5e6,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[1500]*0.6,P_dot_B5e6[1500]*0.75,r'$B =5\times10^{6}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e7,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2000]*0.6,P_dot_B1e7[2000]*0.95,r'$B= 10^{7}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e8,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2300]*0.6,P_dot_B1e8[2300]*0.95,r'$B= 10^{8}$T',fontsize=7,color='black',rotation=-20,zorder=3)
ax_main.plot(P_line,P_dot_B1e9,linestyle='dotted',c='black',zorder=0)
ax_main.text(P_line[2500]*0.6,P_dot_B1e9[2500]*0.95,r'$B= 10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
#plt.plot(P_line,P_dot_B4_4e9,linestyle='dotted',c='black',zorder=0)
#plt.text(P_line[2500]*0.6,P_dot_B4_4e9[2500]*0.95,r'$B= 4.4\times10^{9}$T',fontsize=7,color='black',rotation=-20,zorder=3)
plt.savefig(f'inference_plots/contour_map_ppdot_xray.pdf',dpi=300)

#T histogram
plt.figure(8)
plt.hist(np.log10(all_T),bins=20,edgecolor='red',color='red',histtype='step',alpha=0.5,density=True,label=r'Simulation')
plt.legend()
#plt.xscale('log')
plt.xlabel(r'Temperature of the hot spot in logscale (K)')
plt.ylabel('p.d.f')
plt.savefig('inference_plots/histo_T.pdf',dpi=300)
plt.close()

#r_h histogram
plt.figure(9)
plt.hist(all_rh,bins=20,edgecolor='red',color='red',histtype='step',alpha=0.5,density=True,label=r'Simulation')
plt.legend()
plt.xlabel(r'Radius of the hot spot (m)')
plt.ylabel('p.d.f')
plt.savefig('inference_plots/histo_rh.pdf',dpi=300)
plt.close()

#Log(P) histogram
plt.figure(10)
plt.hist(np.log10(P_sim_all),bins=20,range=(-2,1.5),edgecolor='red',color='red',alpha=1,density=True,label='Simulation',histtype='step')
plt.hist(np.log10(P_obs_all),bins=20,range=(-2,1.5),edgecolor='blue',color='blue',alpha=1,density=True,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel(r'Log($P$) ($P$ in s)')
plt.ylabel('p.d.f')
plt.savefig('inference_plots/histo_period.pdf',dpi=300)
plt.close()

#Log(P_dot) histogram
plt.figure(11)
plt.hist(np.log10(P_dot_sim_all),bins=20,range=(-20,-10),edgecolor='red',color='red',alpha=1,density=True,label='Simulation',histtype='step')
plt.hist(np.log10(P_dot_obs_all),bins=20,range=(-20,-10),edgecolor='blue',color='blue',alpha=1,density=True,label='ATNF data',histtype='step') #Canonical pop
plt.legend()
plt.xlabel(r'Log ($\dot{P}$) ($\ \dot{P}$ in $s.s^{-1})$')
plt.ylabel('p.d.f')
plt.savefig('inference_plots/histo_pdot.pdf',dpi=300)
plt.close()

#cos(alpha0) and cos(alpha) histogram
plt.figure(12)
plt.hist(cos_alpha_all,bins=20,range=(-1,1),edgecolor='red',color='red',alpha=0.5,density=True,histtype='step',label=r'cos($\chi$)')
plt.hist(cos_alpha0_all,bins=20,range=(-1,1),edgecolor='blue',color='blue',alpha=0.5,density=True,histtype='step',label=r'cos($\chi_0$)')
plt.legend()
plt.xlabel(r'cos($\chi$) and cos($\chi_0$)')
plt.ylabel('p.d.f')
plt.savefig('inference_plots/histo_cosalpha.pdf')
plt.close()
