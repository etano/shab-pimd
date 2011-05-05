import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import sys

# Input File
inputFile = []
for i in range(0, len(sys.argv)-1):  
  inputFileName = "inputs/" + str(sys.argv[i+1])
  inputFile.append([])
  inputFile[i] = open(inputFileName,'r')

# Input Data
inputFileLabel = []
scalarData = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  scalarData.append([])
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
      scalarFilePath = "data/traces/scalarTrace-" + fileExtension + ".dat"    
      print "\nReading data from " + scalarFilePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(scalarFilePath)
      scalarData[i].append(CalcStatistics.getAndOutputStats(myArray, myArrayHeadings))
      Plotting.makePlots(myArray, myArrayHeadings, fileExtension)



# THIS IS HARD CODED
# NEEDS TO BE CHANGED FOR X-VALUES
# OTHER THAN 1 - 10

# Beta Values
beta = []
for i in range(0, len(scalarData)): 
  # For every input file 
  beta.append([])
  for j in range(0, len(scalarData[i])):
    # For every beta point
    beta[i].append(j+1)



# Rotate Data into Columns
col = []
for i in range(0, len(scalarData)):
  # For every input file
  col.append([])
  for j in range(0, len(scalarData[i][0])):
    # For every observable
    col[i].append([])
    for k in range(0, len(scalarData[i][0][0])):
      # For every statistic
      col[i][j].append([])
      for l in range(0, len(scalarData[i])):
        # For every beta point
        col[i][j][k].append([])
        col[i][j][k][l] = scalarData[i][l][j][k]

# Generate Plots
print "\nGenerating Plots:"
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.errorbar(beta[j], col[j][i][0], col[j][i][3], label=inputFileLabel[j])
  leg = plt.legend(loc='best')
  # matplotlib.text.Text instances
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle(myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/" + myArrayHeadings[i+1] + "vBeta" + fileExtension + ".png")
  plt.clf()
  print "\nPlot data/figures/" + myArrayHeadings[i+1] + "vBeta" + fileExtension + ".png Generated!"

print "\nDone.\n"
