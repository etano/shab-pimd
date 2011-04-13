import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

nPart = 1
nBead = 100
nD = 1
beta = 10
dt = 0.01
L = 1
useStage = 1
useNH = 1
nNH = 3
nSY = 5

filename = "-" + str(nPart) + "-" + str(nD) + "-" + str(nBead) + "-" + str(beta) + "-" + str(dt) + "-" + str(L) + "-" + str(useStage) + "-" + str(useNH) + "-" + str(nNH) + "-" + str(nSY)
    
(myArray, myArrayHeadings) = ReadData.loadAscii("data/traces/scalarTrace" + filename + ".dat")
#CalcStatistics.getAndOutputStats(myArray, myArrayHeadings)
Plotting.makePlots(myArray, myArrayHeadings, filename)

print "\ndone.\n"
