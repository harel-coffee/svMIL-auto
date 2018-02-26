#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	The purpose of this script is to get more insight into which threshold is the best option for reducing the number of Hi-C interactions.
	Thresholds to look at:
		- Uniform threshold. Simply remove all interactions with a count lower than a specific value
		- Percentage-based threshold. Only retain the interactions with an interaction count higher than the top percentage of interaction counts. 
	There is currently a function to actually filter the data from this script, but there is an additional script (thresholdHiCData.py) that is dedicated to this task. 
"""


### Imports ###

import numpy as np
import time
import pickle
import os
import glob

### Code ###

#Determine if the actual filtering file needs to be made (True) or not, or if we want to create the visualization file (False)
filtering = False

#Determine the threshold
threshold = 5
uniform = True
thresholdPercentage = float(2)

#Specify if this intrachromosomal or interchromosomal dat
intrachromosomal = True

#First, read the interaction data
#I will read the interactions from the pkl file generated in the plotInteractionCountDistribution.py file to not have to re-read the data from disk each time

#Show some statistics on how many interactions we retained
pkl_file = open('interactionCounts_intrachromosomal.pkl', 'rb')

interactionCounts = pickle.load(pkl_file)

#Apply a threshold to the interactionCounts

if uniform == True:
	
	filteredInteractionCounts = dict()
	filteredInteractionCountValue = 0
	unfilteredInteractionCountValue = 0
	for interactionCount in interactionCounts:
		
		if interactionCount >= threshold:
			filteredInteractionCounts[interactionCount] = interactionCounts[interactionCount]
			filteredInteractionCountValue += interactionCounts[interactionCount]
		unfilteredInteractionCountValue += interactionCounts[interactionCount]
	#Show the total value count
	
	print "uniform threshold of ", threshold, ": ", filteredInteractionCountValue
	print "unfiltered: ", unfilteredInteractionCountValue
	
	topInteractions = filteredInteractionCounts
	usedThreshold = threshold
else:
	
	#use a percentage threshold
	
	#first rank all interaction counts
	#Obtain the top x%.
	
	
	#Obtain the unique interaction count values. 
	
	interactionCountKeys = np.array(interactionCounts.keys())
	print "total number of unique interaction counts: ", len(interactionCountKeys)
	sortedInteractionCounts = np.sort(interactionCountKeys)
	
	cutoff = int(round(len(interactionCountKeys) * (thresholdPercentage/100)))
	print "cutoff: ", cutoff
	indices = np.arange(cutoff, len(interactionCountKeys), dtype=int)
	
	
	topInteractions = sortedInteractionCounts[indices]
	
	print topInteractions
	
	filteredInteractionCounts = dict()
	filteredInteractionCountValue = 0
	for interactionCount in topInteractions:
		print "interaction: ", interactionCount
		filteredInteractionCounts[interactionCount] = interactionCounts[interactionCount]
		filteredInteractionCountValue += interactionCounts[interactionCount]
	usedThreshold = thresholdPercentage

print topInteractions
#Show the total value count

print "percentage filtering with threshold of :", usedThreshold, ": ", filteredInteractionCountValue

#To validate if the threshold does not remove too much of the off-diagonal elements, we can visualize the resulting matrix. It may also be that we end up with different thresholds between the interchromosomal
#and intrachromosomal data. 

#Get back the original interactions above the lowest frequency count, and keep all these in a file.

lowestRequiredInteractionCount = min(topInteractions)
print "required threshold: ", lowestRequiredInteractionCount

if filtering == True:

	#Read the original data.
	
	#Filter for when the interaction count is lower than our threshold.
	
	#Write the data to a filtered file if it meets the requirements.

	hiCFolder = "../../../../../Data/HiC/HUVEC_intrachromosomal/5kb_resolution_intrachromosomal/"
	
	#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
	filesToRead = []
	for folderName in os.listdir(hiCFolder):
		#Skip if it is an osx specific file and not a folder
		if os.path.isdir(hiCFolder + "/" + folderName):
			filesToRead.append(hiCFolder + "/" + folderName + "/MAPQG0/*.RAWobserved")
	
	startTime = time.time()
	
	interactionCounts = dict()
	
	fileCount = 0
	for hiCFilePath in filesToRead:
		
		hiCFile = glob.glob(hiCFilePath) #get the actual file path instead of the wildcard one
		print hiCFile[0] #There will always be only one file that matches
		
		currentEdges = [] #refresh the list to make sure that not everything is in memory
		with open(hiCFile[0], "r") as f:
			
			splitFilePath = hiCFile[0].split("/")
			fileName = splitFilePath[len(splitFilePath)-1]
			
			outFile = "thresholdedHiC/" + fileName + "_thresholded_" + str(usedThreshold)
			
			#extract the chromosome name from the file (this is very quick and dirty but works for now)
			splitFileName = fileName.split(".")
			splitFileName = splitFileName[0]
			splitFileName = splitFileName.split("_")
			
			chr1 = splitFileName[0]
			chr2 = 'chr' + str(splitFileName[1])
			
			if intrachromosomal == True:
				chr2 = chr1

			with open(outFile, 'w') as outF:
				
				for line in f:
					line = line.strip()
					splitLine = line.split("\t")
				
					count = int(float(splitLine[2]))
					if count >= lowestRequiredInteractionCount:
						
						if filtering == True: #This is for filtering
							outF.write(line)
							outF.write("\n")
							
						else: #This is for generating an input format that juicer can use to convert it to .hiC format, which can be read by the juicebox visualizer
							
							for interaction in range(0, count): #Do this multiple times, because the hiC format thinks each line is a read, so how often we count an interaction, while the RAWobserved is already grouped
								#Write line to new file
								#Format that we need to convert the data to:
								#str1	chr1	pos1	frag1	str2	chr2	pos2	frag2
								str1 = "0" #this is not used by the visualization toolbox, so we can dummy it
								str2 = str1
								frag1 = str(0)
								frag2 = str(1)
								pos1 = str(splitLine[0])
								pos2 = str(splitLine[1])
								newFileLine = str1 + " " + chr1 + " " + pos1 + " " + frag1 + " " + str2 + " " + chr2 + " " + pos2 + " " + frag2 + "\n"
								outF.write(newFileLine)
						
					
					
						
						
	
		fileCount += 1
		
		

