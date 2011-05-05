#!/bin/sh
varCount="1"
vars=""
for var in "$@"
do
  if [ $varCount -eq "1" ]
  then
    nBins=$var
    varCount=2
  elif [ $varCount -eq "2" ]
  then
    label=$var
    varCount=3
  else
    vars+="$var "
  fi
done
python ScalarPlots.py $vars
python DensityPlots.py $nBins $label $vars
python GrrPlots.py $nBins $label $vars
