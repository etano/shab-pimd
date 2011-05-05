import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import * 
from histogram import *
import sys

# Number of Bins
nBins = str(sys.argv[1])

# Output Figure Label
outputLabel = str(sys.argv[2])

# Input File
inputFile = []
for i in range(0, len(sys.argv)-3):  
  inputFileName = "inputs/" + str(sys.argv[i+3])
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
    plt.hist(RDenData[i][j], int(nBins), normed=True, histtype='step', label=inputFileLabel[i])
leg = plt.legend(loc='best')
# matplotlib.text.Text instances
for t in leg.get_texts():
    t.set_fontsize('xx-small')    # the legend text fontsize
plt.suptitle("Density", fontsize=12)
plt.savefig("data/figures/RDen-" + outputLabel + ".png")
plt.clf()
print "\nPlot data/figures/RDen-" + outputLabel + ".png Generated!"

print "\nDone.\n"
