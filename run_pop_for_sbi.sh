#!/bin/sh
#param1=$1
#param2=$2
#param3=$3
#param4=$4
#param5=$5
#param6=$6
#param7=$7
#param8=$8
#param9=$9
#param10=${10}
#param11=${11}
#param12=${12}
#param13=${13}
#param14=${14}
echo "Run of the normal pulsar population synthesis in progress" 
./YoungPop "$@"
sed -i 's/nan/0.000000e+00/g' wint.txt
python3 plots.py
python3 x_ray_analysis.py
echo "The simulation was run successfully"
