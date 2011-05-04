import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import sys

# Output Figure Label
outputLabel = str(sys.argv[1])

# Input File
inputFile = []
for i in range(0, len(sys.argv)-2):  
  inputFileName = "inputs/" + str(sys.argv[i+2])
  inputFile.append([])
  inputFile[i] = open(inputFileName,'r')

# Input Data
inputFileLabel = []
GrrData = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  GrrData.append([])
  for line in inputFile[i]:
    if firstLine:
    # First line label
      inputFileLabel[i] = str(line)
      firstLine = False
    else:
    # For every beta point
      inputLine = str(line)
      fileExtension = inputLine.replace(' ','-')
      fileExtension = fileExtension.replace('\n','')
      GrrFilePath = "data/traces/GrrTrace-" + fileExtension + ".dat"    
      print "\nReading data from " + GrrFilePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(GrrFilePath)
      GrrData[i].append(myArray)
      #Plotting.makePlots(myArray, myArrayHeadings, fileExtension)

# Generate Plots
print "\nGenerating Plots:"
plt.xlabel("r")
plt.ylabel("Grr")
for i in range(0, len(inputFile)):
  # For every input file 
  for j in range(0, len(GrrData[i])):
    # For every line in the input file
    plt.hist(GrrData[i][j], 100, label=inputFileLabel[i])
plt.legend()
plt.suptitle("Pair Correlation, Grr", fontsize=12)
plt.savefig("data/figures/Grr-" + fileExtension + ".png")
plt.clf()
print "\nPlot data/figures/Grr-" + fileExtension + ".png Generated!"

print "\nDone.\n"
