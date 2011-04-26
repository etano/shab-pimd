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
stats = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  stats.append([])
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
      filePath = "data/traces/scalarTrace-" + fileExtension + ".dat"    
      print "\nReading data from " + filePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(filePath)
      stats[i].append(CalcStatistics.getAndOutputStats(myArray, myArrayHeadings))
      #Plotting.makePlots(myArray, myArrayHeadings, fileExtension)

# Beta Values
beta = []
for i in range(0, len(stats)): 
  # For every input file 
  beta.append([])
  for j in range(0, len(stats[i])):
    # For every beta point
    beta[i].append(j+1)

# Rotate Data into Columns
col = []
for i in range(0, len(stats)):
  # For every input file
  col.append([])
  for j in range(0, len(stats[i][0])):
    # For every observable
    col[i].append([])
    for k in range(0, len(stats[i][0][0])):
      # For every statistic
      col[i][j].append([])
      for l in range(0, len(stats[i])):
        # For every beta point
        col[i][j][k].append([])
        col[i][j][k][l] = stats[i][l][j][k]

# Generate Plots
print "\nGenerating Plots:"
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.errorbar(beta[j], col[j][i][0], col[j][i][3], label=inputFileLabel[j])
  plt.legend()
  plt.suptitle(myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/" + outputLabel + "-" + myArrayHeadings[i+1] + "VsBeta.png")
  plt.clf()
  print "\nPlot data/figures/" + outputLabel + "-" + myArrayHeadings[i+1] + "VsBeta.png Generated!"

print "\nDone.\n"
