#!/bin/sh
start_time=$(date +%s)
#nvcc -g -O0 -fmad=false -rdc=true --ptxas-options=-v main.cu initialize.cu birth_pulsars.cu detection.cu ism_scattering.cu galac_pot.cu evolution.cu ymw16.a -o  YoungPop -lm -lgsl -lgslcblas -lc -lcudadevrt
nvcc -g -O0 -rdc=true -fmad=false -Xcompiler -fno-fast-math --ptxas-options=-v main.cu initialize.cu birth_pulsars.cu detection.cu galac_pot.cu evolution.cu ymw16.a -o YoungPop -lm -lgsl -lgslcblas -lc -lcudadevrt
./YoungPop
sed -i 's/nan/0.000000e+00/g' wint.txt
python3 plots.py
python3 x_ray_analysis.py
end_time=$(date +%s)
elapsed_time=$(( end_time - start_time ))
hours=$(( elapsed_time / 3600 ))
minutes=$(( (elapsed_time % 3600) / 60 ))
seconds=$(( elapsed_time % 60 ))
echo "Elapsed time : $(printf "%02d" $hours):$(printf "%02d" $minutes):$(printf "%02d" $seconds)"
