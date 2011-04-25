import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import sys

inputFileName = str(sys.argv[1])
inputFile = open(inputFileName,'r')

for line in inputFile: 
  inputLine = str(line)
  fileExtension = inputLine.replace(' ','-')
  fileExtension = fileExtension.replace('\n','')
  filePath = "data/traces/scalarTrace-" + fileExtension + ".dat"    
  print "\nReading data from " + filePath + "."
  (myArray, myArrayHeadings) = ReadData.loadAscii(filePath)
  #CalcStatistics.getAndOutputStats(myArray, myArrayHeadings)
  Plotting.makePlots(myArray, myArrayHeadings, fileExtension)

print "\nDone.\n"
