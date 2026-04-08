import healpy as hp
import numpy as np
import re
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC

#Remove old data
os.system(f'rm temp.txt')

# Load the coordinates data
reg_4=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
with open("l_b_coord_sim.txt","r") as f:
    data_lb=re.findall(reg_4,f.read())

number_of_sources=int(len(data_lb)/2)

l,b=[],[]
for i in range(number_of_sources):
    l.append(float(data_lb[2*i]))
    b.append(float(data_lb[2*i+1]))

#Example test
#l=[89,21,4]
#b=[2,-39,72]
#number_of_sources=3
# Load healpix map
sky_map = hp.read_map('haslam408_dsds_Remazeilles2014_ns2048.fits')  # Remplacez par le nom du fichier téléchargé

#Get temperature of the sky
nside = hp.get_nside(sky_map)
theta = np.radians(90 - np.array(b))
phi = np.radians(np.array(l))       
pix = hp.ang2pix(nside, theta, phi)
temperatures = sky_map[pix]
np.savetxt('temp.txt', temperatures, fmt='%.6e')

# Convert equatorial coordinates in galactic coordinates
#for i in range(number_of_sources):
    # Get the corresponding pixel
    #nside = hp.get_nside(sky_map)
    #theta = np.radians(90 - b[i])
    #phi = np.radians(l[i])
    #pix = hp.ang2pix(nside, theta, phi)
    # Extract sky temperature
    #temperature = sky_map[pix]
    #with open('temp.txt','a') as f:
    #    f.write(f'{temperature}\n')
    #Test if this region was observed by either XMM-Newton or Chandra
    #X_obs_XMM=int(xmm_observed(l[i],b[i]))
    #X_obs_chandra=int(chandra_observed(l[i],b[i]))
    #with open('sky_X_obs_XMM.txt','a') as f:
    #    f.write(f'{X_obs_XMM}\n')
    #with open('sky_X_obs_chandra.txt','a') as f:
    #    f.write(f'{X_obs_chandra}\n')
#print(f"Température du ciel à 408 MHz : {temperature:.2f} K")
