#!/bin/sh
nvcc -g main.cu initialize.cu birth_pulsars.cu detection.cu ism_scattering.cu galac_pot.cu evolution.cu  -o  YoungPop -lm -lgsl -lgslcblas  


##gcc -g -W -Wall -L/usr/local/gsl-2.5 young_pulsar_pop.c initialize.c birth_pulsars.c evolution.c -o YoungPop -lm -lgsl -lgslcblas -lchealpix

