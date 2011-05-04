import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import sys

# Number of Bins
nBins = str(sys.argv[1])

# Input File
inputFile = []
for i in range(0, len(sys.argv)-2):  
  inputFileName = "inputs/" + str(sys.argv[i+2])
  inputFile.append([])
  inputFile[i] = open(inputFileName,'r')

# Input Data
inputFileLabel = []
RDenData = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  RDenData.append([])
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
      RDenFilePath = "data/traces/RDenTrace-" + fileExtension + ".dat"    
      print "\nReading data from " + RDenFilePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(RDenFilePath)
      RDenData[i].append(myArray)

# Generate Plots
print "\nGenerating Plots:"
plt.xlabel("r")
plt.ylabel("Density")
for i in range(0, len(inputFile)):
  # For every input file 
  for j in range(0, len(RDenData[i])):
    # For every line in the input file
    plt.hist(RDenData[i][j], int(nBins), normed=True, label=inputFileLabel[i])
plt.legend()
plt.suptitle("Density", fontsize=12)
plt.savefig("data/figures/RDen-" + fileExtension + ".png")
plt.clf()
print "\nPlot data/figures/RDen-" + fileExtension + ".png Generated!"

print "\nDone.\n"
