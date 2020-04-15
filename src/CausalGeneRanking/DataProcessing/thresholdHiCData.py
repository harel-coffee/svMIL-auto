#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	The purpose of this script is to filter the HiC data given a defined threshold. Every interaction with a measured interaction count above the defined threshold will be written
	to the defined new file. 
"""


### Imports ###

import numpy as np
import matplotlib.pyplot as plt
import os
import re
import time
import glob
import pickle
import sys
### Settings ###

hiCFolder = sys.argv[1] #Main folder from which all directories will be read
fileExtension = "/MAPQG0/*.RAWobserved" #In the directories of the folder specified above, we further search at this path for files with the fiven extension
interactionType = "intrachromosomal"

outDir = sys.argv[2]

finalOutDir = outDir + "/" + interactionType + "/"
print(finalOutDir)

threshold = 5 #For now I only use a uniform threshold, other thresholds could be implemented

pickleExists = False #use this switch to use a pickle file from a previous run if already generated, then the data does not need to be read again


### Code ###


### 1. Input: read all interactions from a provided location into a dictionary

if pickleExists == True:

	pkl_file = open('interactionCounts_' + interactionType + '.pkl', 'rb') 
	interactionCounts = pickle.load(pkl_file)

else:
	
	#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
	filesToRead = []
	for folderName in os.listdir(hiCFolder):
		#Skip if it is an osx specific file and not a folder
		if os.path.isdir(hiCFolder + "/" + folderName):
			filesToRead.append(hiCFolder + "/" + folderName + fileExtension)
	
	startTime = time.time()
	
	interactionCounts = dict()

	fileCount = 0
	for hiCFilePath in filesToRead:
		hiCFile = glob.glob(hiCFilePath) #get the actual file path instead of the wildcard one
		print(hiCFile[0]) #There will always be only one file that matches
		
		currentEdges = [] #refresh the list to make sure that not everything is in memory
		with open(hiCFile[0], "r") as f:
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")

				count = int(float(splitLine[2]))
				if count > 0:
					
					if count in interactionCounts:
						interactionCounts[count] += 1
					else:
						interactionCounts[count] = 1
		
		fileCount += 1
		
	endTime = time.time()
	print("Took ", endTime - startTime, " seconds to store edges from all files")
	
	output = open('interactionCounts_' + interactionType + '.pkl', 'wb')
	
	# Pickle dictionary using protocol 0.
	pickle.dump(interactionCounts, output)


### 2. Do the filtering given our defined threshold
filteredInteractionCounts = dict()
filteredInteractionCountValue = 0
unfilteredInteractionCountValue = 0
for interactionCount in interactionCounts:
	
	if interactionCount >= threshold:
		filteredInteractionCounts[interactionCount] = interactionCounts[interactionCount]
		filteredInteractionCountValue += interactionCounts[interactionCount]
	unfilteredInteractionCountValue += interactionCounts[interactionCount]
#Show the total value count

topInteractions = filteredInteractionCounts

### 3. Write the interactions to the provided output location and filter by the defined threshold


#Prepare the output directory
if not os.path.exists(finalOutDir):
   os.makedirs(finalOutDir)


lowestRequiredInteractionCount = min(topInteractions)
print("required threshold: ", lowestRequiredInteractionCount)

#Read the original data. Here we do need to actually read the files again. 

#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
filesToRead = []
for folderName in os.listdir(hiCFolder):
	#Skip if it is an osx specific file and not a folder
	if os.path.isdir(hiCFolder + "/" + folderName):
		filesToRead.append(hiCFolder + "/" + folderName + fileExtension)

startTime = time.time()

interactionCounts = dict()

fileCount = 0
for hiCFilePath in filesToRead:
	
	hiCFile = glob.glob(hiCFilePath) #get the actual file path instead of the wildcard one
	print(hiCFile[0]) #There will always be only one file that matches
	
	currentEdges = [] #refresh the list to make sure that not everything is in memory
	with open(hiCFile[0], "r") as f:
		
		splitFilePath = hiCFile[0].split("/")
		fileName = splitFilePath[len(splitFilePath)-1]
		
		outFile = finalOutDir + '/' + fileName + "_thresholded_" + str(threshold)
		print(outFile)
		
		#extract the chromosome name from the file (this is very quick and dirty but works for now)
		splitFileName = fileName.split(".")
		splitFileName = splitFileName[0]
		splitFileName = splitFileName.split("_")
		
		chr1 = splitFileName[0]
		chr2 = 'chr' + str(splitFileName[1])
		
		if interactionType == "intrachromosomal":
			chr2 = chr1

		with open(outFile, 'w') as outF:
			
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
			
				count = int(float(splitLine[2]))
				if count >= lowestRequiredInteractionCount:

					outF.write(line)
					outF.write("\n")
						

	fileCount += 1
		
		

