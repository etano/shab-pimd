#!/bin/sh
python ScalarPlots.py $1
python DensityPlots.py $2 $1
python GrrPlots.py $2 $1
