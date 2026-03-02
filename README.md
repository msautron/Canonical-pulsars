# User guide

This code generates a synthetic population of normal pulsars, considering radio and gamma-ray emission. If you download this repository, you are comparing the simulation to the detection made by Parkes, FAST and Fermi/LAT. 
A GPU is needed. The library GSL is needed.

## Compilation 

Use bash run_pop.sh to compile, execute the code and generate the plots ! 

## Description of the code

`macro.h` -> Contains the declaration of most variables and lists.

`initialize.cu` -> Declare also some variables + allows to finalize the declaration of the tables. The parameters can be changed here. 

`main.cu` -> Spatial evolution, made with GPU is here, allows to choose the order of use of each function. Free the memory at the end.

`cn.h`, `dmdtau.cu`, `dora.cu`, `fermibubble.cu`, `frb_d.cu`, `galcen.cu`, `gum.cu`, `lmc.cu`, `localbubble.cu`, `ne_crd.cu`, `nps.cu`, `smc.cu`, `spiral.cu`, `spiral.txt`, `thick.cu`, `thin.cu`, `ymw16`, `ymw16a`, `ymw16.cu`, `ymw16_ne`, `ymw16_ne.cu`, `ymw16par.cu` and `ymw16par.txt` -> files related to ISM dispersion of the width of the radio pulse profile. Adapted to GPU from Yao et al. (2017)

`galac_pot.cu` -> Contains function for initial speed.

`3PC_Catalog_20230803.xls` and `3PC_SensitivityMap_20230629.fits` -> Contains the 3PC catalog and Fermi/LAT sensitivity map

`birth_pulsars.cu` -> Contains the initial distribution for P, B and the age.

`evolution.cu` -> Contains the functions to evolve Omega, Omega_dot and B. 

`detection.cu` -> Handles the detection in radio + gamma. Most important functions are there.

`facteur_omega.dat` -> File obtained from J.Pétri allowing to obtain the anisotropy factor from simulation (cf J.Pétri (2024)). 

`get_temp.py` -> Python file allowing to compute the temperature of the sky at a position in the sky. 

`plots.py` -> Makes the plots. Look here if you want to know which file is needed to obtain plots.

`sensitivity_3PC.py` -> Compute sensitivity in gamma thanks to the fits map. 

`fast_fermi_pmps.txt` -> Contains the data to compare with normal pulsars detected by FAST, Fermi and Parkes. 
