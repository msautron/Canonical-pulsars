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
num_sim=200
num_training=10
num_check=5

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
#Repeat_sim_and_save(num_sim,num_training,prior)
os.system(f'rm P_Pdot_positions.txt')
os.system(f'rm PPdot_alpha_all.txt')
os.system(f'rm PF_check.txt')
os.system(f'rm x_file.txt')
os.system(f'rm x_file2.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm gamma_peak_sep.txt')
os.system(f'rm Fg_flux.txt')
os.system(f'rm xi_rho_data.txt')
os.system(f'rm nb_orbit.txt')
os.system(f'rm wint.txt')
os.system(f'rm wr.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm init_P_B.txt')
os.system(f'rm Lgamma.txt')
posterior,observation=SBI_from_datafiles(prior,valid_or_obs,show)

#Using the parameters of SBI for num_check simulations
samples = posterior.sample((num_check,),x=observation)
os.system(f'rm P_Pdot_positions.txt')
os.system(f'rm PPdot_alpha_all.txt')
os.system(f'rm PF_check.txt')
os.system(f'rm x_file.txt')
os.system(f'rm x_file2.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm gamma_peak_sep.txt')
os.system(f'rm Fg_flux.txt')
os.system(f'rm xi_rho_data.txt')
os.system(f'rm nb_orbit.txt')
os.system(f'rm wint.txt')
os.system(f'rm wr.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm init_P_B.txt')
os.system(f'rm Lgamma.txt')
params = {
            'sigma_b'    : float(samples[0,0]),
            'b_mean'     : float(10**(samples[0,1])),
            'p_mean'     : float(samples[0,2]*1e-3),
            'sigma_p'    : float(samples[0,3]),
            'BR'    : float(np.round(samples[0,4])),
            'A_propto1' : float(10**(samples[0,5])),
            'D_propto1' : float(10**(samples[0,6])),
            'M_for_K' : float(samples[0,7]),
            'R_for_K' : float(samples[0,8]),
            'tau_d'   : float((10**samples[0,9])*365*24*3600),
            'alpha_d' : float(samples[0,10])
            }
print(params)
cmd = [
    'bash',
    'run_pop_post_inference.sh',
    str(params['sigma_b']),
    str(params['b_mean']),
    str(params['p_mean']),
    str(params['sigma_p']),
    str(params['BR']),
    str(params['A_propto1']),
    str(params['D_propto1']),
    str(params['M_for_K']),
    str(params['R_for_K']),
    str(params['tau_d']),
    str(params['alpha_d'])
]
print("Command executed: ")
print(" ".join(cmd))
sim_run=subprocess.run([f'bash','run_pop_post_inference.sh', str(params['sigma_b']), str(params['b_mean']), str(params['p_mean']), str(params['sigma_p']), str(params['BR']), str(params['A_propto1']), str(params['D_propto1']), str(params['M_for_K']), str(params['R_for_K']), str(params['tau_d']), str(params['alpha_d'])], capture_output=True,text=True)
print("STDOUT :", sim_run.stdout)
print("STDERR :", sim_run.stderr)
print("Return code :", sim_run.returncode)
#os.system(f"sed -i 's/nan/0.000000e+00/g' wint.txt")
os.system(f"python3 plots.py")
os.system(f'rm P_Pdot_positions.txt')
os.system(f'rm PPdot_alpha_all.txt')
os.system(f'rm PF_check.txt')
os.system(f'rm x_file.txt')
os.system(f'rm x_file2.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm gamma_peak_sep.txt')
os.system(f'rm Fg_flux.txt')
os.system(f'rm xi_rho_data.txt')
os.system(f'rm nb_orbit.txt')
os.system(f'rm wint.txt')
os.system(f'rm wr.txt')
os.system(f'rm info_supp.txt')
os.system(f'rm init_P_B.txt')
os.system(f'rm Lgamma.txt')
for i in range(1,num_check):
    params = {
                'sigma_b'    : float(samples[i,0]),
                'b_mean'     : float(10**(samples[i,1])),
                'p_mean'     : float(samples[i,2]*1e-3),
                'sigma_p'    : float(samples[i,3]),
                'BR'    : float(np.round(samples[i,4])),
                'A_propto1' : float(10**(samples[i,5])),
                'D_propto1' : float(10**(samples[i,6])),
                'M_for_K' : float(samples[i,7]),
                'R_for_K' : float(samples[i,8]),
                'tau_d'   : float((10**samples[i,9])*365*24*3600),
                'alpha_d' : float(samples[i,10])
                }
    sim_run=subprocess.run([f'bash','run_pop_post_inference.sh', str(params['sigma_b']), str(params['b_mean']), str(params['p_mean']), str(params['sigma_p']), str(params['BR']), str(params['A_propto1']), str(params['D_propto1']), str(params['M_for_K']), str(params['R_for_K']), str(params['tau_d']), str(params['alpha_d'])], capture_output=True,text=True)
    print("STDOUT :", sim_run.stdout)
    print("STDERR :", sim_run.stderr)
    print("Return code :", sim_run.returncode)
    #os.system(f"sed -i 's/nan/0.000000e+00/g' wint.txt")
    #os.system(f"python3 plots_post_inference.py")
    print(f"This was the simulation number {i+1}")

os.system(f"python3 plots_post_inference.py")

#Compute the time of execution
end_time=time.time()
elapsed_time = end_time -start_time

hours, remainder = divmod(elapsed_time,3600)
minutes, seconds = divmod(remainder,60)

print(f"It took {int(hours)}:{int(minutes)}:{int(seconds)} to run \n")
