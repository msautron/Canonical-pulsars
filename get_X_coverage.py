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
from astropy.coordinates import match_coordinates_sky

#Remove old data
os.system(f'rm sky_X_obs_XMM.txt')
os.system(f'rm sky_X_obs_chandra.txt')

# Load the coordinates data
reg_4=re.compile("-*\d{1}[.]\d{6}[eE]*[+-]*\d{2}")
with open("l_b_coord_sim.txt","r") as f:
    data_lb_X=re.findall(reg_4,f.read())

number_of_sources=int(len(data_lb_X)/2)

l,b=[],[]
for i in range(number_of_sources):
    l.append(float(data_lb_X[2*i]))
    b.append(float(data_lb_X[2*i+1]))

#Getting info about the sky coverage
#XMM-Newton sky coverage
data_XMM = fits.open("4xmmdr14_240411.fits")[1].data
ra = data_XMM["RA"]
dec = data_XMM["DEC"]
sources = SkyCoord(ra*u.deg, dec*u.deg)

#def xmm_observed(l,b): #One by one
#    pulsar = SkyCoord(l*u.deg, b*u.deg, frame="galactic").icrs
#    sep = pulsar.separation(sources)
#    return np.min(sep) < 30*u.arcmin

#def xmm_observed(l_list, b_list):
#    pulsars = SkyCoord(l=l_list*u.deg, b=b_list*u.deg, frame="galactic").icrs
#    sep_matrix = pulsars[:, None].separation(sources[None, :])
#    min_sep = sep_matrix.min(axis=1)
#    return min_sep < 30*u.arcmin

#def xmm_observed(l_list, b_list, sources, max_sep=30*u.arcmin):
#    results = []
#    for l, b in zip(l_list, b_list):
#        pulsar = SkyCoord(l=l*u.deg, b=b*u.deg, frame="galactic").icrs
#        sep = pulsar.separation(sources)
#        results.append(np.min(sep) < max_sep)
#    return results

def xmm_observed(l_list, b_list, sources, max_sep=30*u.arcmin):
    pulsars = SkyCoord(l=l_list*u.deg, b=b_list*u.deg, frame="galactic").icrs
    idx, sep2d, _ = match_coordinates_sky(pulsars, sources)
    return sep2d < max_sep

X_obs_XMM=xmm_observed(l,b,sources,max_sep=30*u.arcmin)
np.savetxt("sky_X_obs_XMM.txt",X_obs_XMM,fmt="%d")

#Chandra sky coverage
moc_chandra=MOC.from_fits("ChandraMOC11_nograting.fits")

def chandra_observed(l,b):
    pulsar = SkyCoord(l*u.deg, b*u.deg, frame="galactic").icrs
    return moc_chandra.contains(pulsar.ra,pulsar.dec)

X_obs_chandra=chandra_observed(l,b)
np.savetxt("sky_X_obs_chandra.txt",X_obs_chandra,fmt="%d")
