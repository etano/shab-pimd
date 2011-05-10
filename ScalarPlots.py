import ReadData 
import Plotting
import CalcStatistics
import math
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
scalarData = []
beta = []
effLabels = []
effX = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  scalarData.append([])
  beta.append([])
  effLabels.append([])
  effX.append([])
  j = 0
  for line in inputFile[i]:
    if firstLine:
    # First line label
      inputFileLabel[i] = str(line)
      firstLine = False
    else:
    # For every beta point
      inputLine = str(line)

      entry = ''
      params = []
      for k in range(len(line)):
        if line[k]!=' ':
          entry += line[k]
        else:
          params += [entry]
          entry=''
      params += [entry]
      entry=''
      betaLabel = params[3]
      beta[i].append(float(params[3]))
      effX[i].append(j)
      nPart = float(params[0])
      nD = float(params[1])
      interaction = float(params[11])
      effLabels[i].append(params[12] + "," + params[13] + "," + params[14])
      fileExtension = inputLine.replace(' ','-')
      fileExtension = fileExtension.replace('\n','')
      scalarFilePath = "data/traces/scalarTrace-" + fileExtension + ".dat"    
      print "\nReading data from " + scalarFilePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(scalarFilePath)
      scalarData[i].append(CalcStatistics.getAndOutputStats(myArray, myArrayHeadings))
      #Plotting.makePlots(myArray, myArrayHeadings, fileExtension)
      
      j += 1

print "\nSorting Data..."
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

# Normal Plots
print "\nGenerating Plots..."
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.errorbar(beta[j], col[j][i][0], col[j][i][3], label=inputFileLabel[j])

  # Exact Values
  x = arange(1.0,10.0,0.01)
  if (interaction==0):  
    if ((myArrayHeadings[i+1]=="PE") or (myArrayHeadings[i+1]=="VE")):  
      y =  nPart*nD*0.5/tanh(x/2.0)
      yclass = nPart*nD/x
      plt.plot(x, y, label="Exact")
      plt.plot(x, yclass, label="Classical")
    if (myArrayHeadings[i+1]=="R2"):  
      y = nD*0.5/tanh(x/2.0)
      yclass = nD/x
      plt.plot(x, y, label="Exact")
      plt.plot(x, yclass, label="Classical")

  # Legend
  leg = plt.legend(loc='best')
  # matplotlib.text.Text instances
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle(myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

# Scaled Plots
for i in range(0, len(col[0])):
  # For every observable)
  if (interaction==0):  
    if ((myArrayHeadings[i+1]=="PE") or (myArrayHeadings[i+1]=="VE") or (myArrayHeadings[i+1]=="R2")):  
      plt.xlabel("Beta")
      plt.ylabel(myArrayHeadings[i+1])
      for j in range(0, len(col)):
        # For every input file 
        yclass2 = nPart*nD/asarray(beta[j])
        plt.errorbar(beta[j], col[j][i][0]-yclass2, col[j][i][3], label=inputFileLabel[j])

      # Exact Values
      x = arange(min(beta[0]),max(beta[0]),0.01)
      y =  nPart*nD*0.5/tanh(x/2.0)
      yclass = nPart*nD/x
      if (myArrayHeadings[i+1]=="R2"):  
        y = y/nPart
        yclass = yclass/nPart
      plt.plot(x, y-yclass, label="Exact")

      # Legend
      leg = plt.legend(loc='best')
      # matplotlib.text.Text instances
      for t in leg.get_texts():
        t.set_fontsize('xx-small')    # the legend text fontsize
      plt.suptitle(myArrayHeadings[i+1] + " Scaled vs Beta", fontsize=12)
      plt.savefig("data/figures/" + myArrayHeadings[i+1] + "ScaledvBeta-" + outputLabel + ".png")
      plt.clf()
      print "\nPlot data/figures/" + myArrayHeadings[i+1] + "ScaledvBeta-" + outputLabel + ".png Generated!"

# Error Plots
print "\nGenerating Error Plots..."
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.plot(beta[j], col[j][i][3], label=inputFileLabel[j])
  leg = plt.legend(loc='best')
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle("Error in " + myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/Err" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/Err" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

# Efficiency Plots
print "\nGenerating Efficiency Plots..."
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("(Length of N-H Chain, Order of S-Y Expansion, Number of Cycles)")
  plt.ylabel("Effciency in " + myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.bar(effX[j], (1.0/(asarray(col[j][0][0])*asarray(col[j][i][3])*asarray(col[j][i][3]))).tolist(), label=inputFileLabel[j])
  plt.xticks(arange(len(effLabels[0]))+1, effLabels[0], rotation='vertical', fontsize='xx-small')
  plt.suptitle("Efficiency of " + myArrayHeadings[i+1] + " for Beta =" + betaLabel, fontsize=12)
  plt.savefig("data/figures/EffBar" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/EffBar" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

# Efficiency Plots
print "\nGenerating Efficiency Plots..."
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.plot(beta[j], (1.0/(asarray(col[j][0][0])*asarray(col[j][i][3])*asarray(col[j][i][3]))).tolist(), label=inputFileLabel[j])
  leg = plt.legend(loc='best')
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle("Efficiency of " + myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/Eff" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/Eff" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

print "\nDone.\n"
