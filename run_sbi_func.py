import torch
import sbi.utils as utils
import numpy as np
import matplotlib.pyplot as plt
from sbi.inference.base import infer
from sbi.analysis import pairplot
from sbi.inference import prepare_for_sbi
from sbi.inference import simulate_for_sbi
from sbi.inference import SNPE
from sbi_tools_CNN import *
import re
from tensorflow.keras import layers, models
import tensorflow as tf
import subprocess
import os
from astropy.coordinates import SkyCoord
from scipy.stats import qmc,gaussian_kde
from scipy import stats
from scipy.interpolate import interp1d
import time

#Measure time of execution
start_time=time.time()

#Choose number of simulations and training
num_sim=100
num_training=200
#num_check=1000

#Sampling the input parameters with LHS method
sigma_b_prior    = np.array([0,1]) # standard deviation of birth magnetic field
b_mean_prior     = np.log10(np.array([2.0e8,3.0e8])) # mean of birth magnetic field [log (T)]
sigma_p_prior    = np.array([0,1]) # standard deviation of birth period
p_mean_prior     = np.array([1,200]) # mean of birth period [ms] converted in s in simulator
BR_prior         = np.array([25,55]) # BR for the NS
#pcst_prior       = np.array([25,26.5]) # Power of the constant for the Lg law of Kalapotharakos et al. (2019)
#pb_prior         = np.array([0.02,0.2]) # Power associated at B for the Lg law of Kalapotharakos et al. (2019)
#pe_prior         = np.array([0.3,0.7]) # Power associated at Edot for the Lg law of Kalapotharakos et al. (2019)
A_propto_prior1  = np.log10(np.array([1e8,8e8])) # Constant of the relation between T and P,Pdot (Harding & Muslimov (2001))
D_propto_prior1  = np.log10(np.array([1e4,8e4])) # Constant of the relation between r_h and R_NS,R_L (Pétri & Mitra (2020))
M_for_K_prior    = np.array([1.0,2.3]) # Mass of NS (in solar mass)
R_for_K_prior    = np.array([8000,16000]) # Radius of NS (in m)
tau_d_prior      = np.log10(np.array([1e4,1e8])) # Tau_d (log yr)
alpha_d_prior    = np.array([0.1,4]) #alpha_d no units
#D_propto_prior2  = np.log10(np.array([1e4,8e4])) # Constant of the relation between T and P,Pdot (Harding & Muslimov (2001))
#A_propto_prior2  = np.log10(np.array([1e8,8e8])) # Constant of the relation between r_h and R_NS,R_L (Pétri & Mitra (2020))

prior = utils.BoxUniform(
    low  = torch.tensor([sigma_b_prior[0],
                         b_mean_prior[0],
                         p_mean_prior[0],
                         sigma_p_prior[0],
                         BR_prior[0],
                         #pcst_prior[0],
                         #pb_prior[0],
                         #pe_prior[0],
                         A_propto_prior1[0],
                         D_propto_prior1[0],
                         M_for_K_prior[0],
                         R_for_K_prior[0],
                         tau_d_prior[0],
                         alpha_d_prior[0]
                         #D_propto_prior2[0],
                         #A_propto_prior2[0]
                        ]),
    high = torch.tensor([sigma_b_prior[1],
                         b_mean_prior[1],
                         p_mean_prior[1],
                         sigma_p_prior[1],
                         BR_prior[1],
                         #pcst_prior[1],
                         #pb_prior[1],
                         #pe_prior[1],
                         A_propto_prior1[1],
                         D_propto_prior1[1],
                         M_for_K_prior[1],
                         R_for_K_prior[1],
                         tau_d_prior[1],
                         alpha_d_prior[1]
                         #D_propto_prior2[1],
                         #A_propto_prior2[1]
                        ])
)

#Test
#t1,t2,t3=simulator([0.5,np.log10(2.75e8),129,0.45,25,45,80,8000,50000,0.3,0.58,np.log10(2e5),np.log10(5e4),np.log10(7e4),np.log10(2e9),np.log10(1e8),np.log10(3e8),26.15,0.06,0.6,np.log10(1.47e9),np.log10(883.1),1.4,12000])

#Allow to run the inference pipeline with the data stocked in result_inference.txt params_inference.txt params_training.txt result_training.txt
valid_or_obs=True #True -> run validation , False -> run the pipeline with comparison with the observations 
show=False #True -> show the cornerplot, False -> Do not show the corner plot
Repeat_sim_and_save(num_sim,num_training,prior)
#posterior,observation=SBI_from_datafiles(prior,valid_or_obs,show)

#Using the parameters of SBI for num_check simulations
#samples = posterior.sample((num_check,),x=observation)
#for i in range(num_check):
#    params = {
#                'sigma_b'    : float(samples[i,0]),
#                'b_mean'     : float(10**(samples[i,1])),
#                'p_mean'     : float(samples[i,2]*1e-3),
#                'sigma_p'    : float(samples[i,3]),
#                'birth_rate' : float(torch.round(samples[i,6]))
#                }
#    subprocess.run([f'{path_to_data}YoungPop',
#                str(params['sigma_b']), str(params['b_mean']), str(params['p_mean']), str(params['sigma_p']), str(params['alpha_d']), str(params['tau_d']), str(params['birth_rate'])], check=True,cwd=path_to_data)
#    os.system(f"python3 {path_to_data}plots.py {i+1}")
#    print(f"This was the simulation number {i+1}")

#Compute the time of execution
end_time=time.time()
elapsed_time = end_time -start_time

hours, remainder = divmod(elapsed_time,3600)
minutes, seconds = divmod(remainder,60)

print(f"It took {int(hours)}:{int(minutes)}:{int(seconds)} to run \n")
