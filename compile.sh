#!/bin/sh
gcc -g -W -Wall -I/home/dirson/gsl/include  -L/home/dirson/gsl/lib main.c initialize.c birth_pulsars.c evolution.c detection.c ism_scattering.c -o YoungPop -lm -lgsl -lgslcblas  


##gcc -g -W -Wall -L/usr/local/gsl-2.5 young_pulsar_pop.c initialize.c birth_pulsars.c evolution.c -o YoungPop -lm -lgsl -lgslcblas -lchealpix

