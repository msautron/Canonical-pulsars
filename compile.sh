#!/bin/sh
gcc -g -gstabs -W -Wall -I/home/matteo/gsl/include  -L/home/matteo/gsl/lib main.c initialize.c birth_pulsars.c evolution.c detection.c ism_scattering.c galac_pot.c  -o  YoungPop -lm -lgsl -lgslcblas  


##gcc -g -W -Wall -L/usr/local/gsl-2.5 young_pulsar_pop.c initialize.c birth_pulsars.c evolution.c -o YoungPop -lm -lgsl -lgslcblas -lchealpix

