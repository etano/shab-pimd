import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *

nPart = 1
nBead = 400
nD = 1
beta = 15.8
dt = 0.01

filename = "-" + str(nPart) + "-" + str(nD) + "-" + str(nBead) + "-" + str(beta) + "-" + str(dt)
    
(myArray, myArrayHeadings) = ReadData.loadAscii("data/traces/scalarTrace" + filename + ".dat")
#CalcStatistics.getAndOutputStats(myArray, myArrayHeadings)
Plotting.makePlots(myArray, myArrayHeadings, filename)

print "\ndone.\n"
