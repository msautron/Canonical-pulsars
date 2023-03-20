import sys
sys.path.append("/home/matteo/Documents/PFE/galpy")
#from galpy.potential import MiyamotoNagaiPotential
from galpy.potential import MWPotential2014
from galpy.orbit import Orbit
import matplotlib.pyplot as plt
import numpy as np
import re
from astropy import units
import argparse

#Get how many time the pulsar is going to evolve
parser=argparse.ArgumentParser(description='Get the age of the pulsar')
parser.add_argument('-a',type=float,dest='age',action='store')
args=parser.parse_args()

#Constant needed
yr_sec=365*24*60*60
sec_Gyr=(1/yr_sec)*10e-9
Tmilky=int(13.5)
step=int(10e-7)
age=args.age*sec_Gyr
print(args.age)
print(sec_Gyr)
print(age)

#Put the data at the right place
var=''
reg=re.compile("-*.{12}[|]{1}")

with open("n_data_fp.txt","r") as f:
    coord=re.findall(reg,f.read())

for i in range(len(coord)):
    for j in range(len(coord[i])-1):
        var+=coord[i][j]
    coord[i]=float(var)
    var=''

#Integration with the gravitational potential

pot=MWPotential2014
#pot=MiyamotoNagaiPotential(a=0.5,b=0.375,amp=1.,normalize=1.)
o=Orbit([coord[0]*units.kpc,coord[1]*units.km/units.s,coord[2]*units.km/units.s,coord[3]*units.kpc,coord[4]*units.km/units.s,coord[5]*units.rad])
ts=np.linspace(int(Tmilky-age),int(Tmilky),1000000)*units.Gyr
o.integrate(ts,pot,method='symplec6_c')

file=open("data_for_c.txt","a")
file.write(f"{o.x(Tmilky*units.Gyr)} ")
file.write(f"{o.y(Tmilky*units.Gyr)} ")
file.write(f"{o.z(Tmilky*units.Gyr)} ")
file.write(f"{o.vx(Tmilky*units.Gyr)} ")
file.write(f"{o.vy(Tmilky*units.Gyr)} ")
file.write(f"{o.vz(Tmilky*units.Gyr)} \n")
file.close()
