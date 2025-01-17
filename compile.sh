#!/bin/sh
#nvcc -g -O0 -fmad=false -rdc=true --ptxas-options=-v main.cu initialize.cu birth_pulsars.cu detection.cu ism_scattering.cu galac_pot.cu evolution.cu ymw16.a -o  YoungPop -lm -lgsl -lgslcblas -lc -lcudadevrt
nvcc -g -O0 -rdc=true -fmad=false -Xcompiler -fno-fast-math --ptxas-options=-v main.cu initialize.cu birth_pulsars.cu detection.cu ism_scattering.cu galac_pot.cu evolution.cu ymw16.a -o YoungPop -lm -lgsl -lgslcblas -lc -lcudadevrt

