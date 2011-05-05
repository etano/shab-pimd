import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
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
    plt.hist(GrrData[i][j], int(nBins), normed=True, histtype='step', label=inputFileLabel[i])
leg = plt.legend(loc='best')
# matplotlib.text.Text instances
for t in leg.get_texts():
    t.set_fontsize('xx-small')    # the legend text fontsize
plt.suptitle("Pair Correlation, Grr", fontsize=12)
plt.savefig("data/figures/Grr-" + outputLabel + ".png")
plt.clf()
print "\nPlot data/figures/Grr-" + outputLabel + ".png Generated!"

print "\nDone.\n"
