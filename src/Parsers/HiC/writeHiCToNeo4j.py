#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This script can read HiC matrices (both interchromosomal and intrachromosomal) in the format of the HUVEC HiC matrices from the Rao et al 2014 paper (see GEO GSE63525). 
	
	If you want to run the script for the intrachromosomal data, you will have to fix the paths yourself for now. This will be fixed later. 
	
	Important:
		- The HiC data is not stored on GitHub as these files are large. Make sure to download these to the folder specified below. These paths will need to be made non-hardcoded later to overcome this issue. 
	
	TODO:
		- Improve documentation of the code
		- Remove hardcoded paths

"""

### Imports ###

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import glob
import pickle
import re
from subprocess import call

### Code ###

#Read the edges and directly write these to a csv file that can be read by Neo4J. Using Networx or other libraries runs into memory issues. 

hiCFolder = "../../../../../Data/HiC/HUVEC_interchromosomal/5kb_resolution_interchromosomal/" #Replace this path if you want to run the script yourself

#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
filesToRead = []
for folderName in os.listdir(hiCFolder):
	#Skip if it is an osx specific file and not a folder
	if os.path.isdir(hiCFolder + "/" + folderName):
		filesToRead.append(hiCFolder + "/" + folderName + "/MAPQG0/*.RAWobserved")

startTime = time.time()

#The output files to write the results to, these will be input for Neo4J
regionsFile = "regions.csv"
relationsFile = "regions_regions_rel.csv"

subdir = "Interchromosomal/"

if not os.path.exists(subdir):
   os.makedirs(subdir)


seenPositions = dict() #Keep a dictionary with previous regions that were already seen. Neo4J has issues with resolving collisions when there are many non-unique regions, so filtering here saves computational time. 

previousLen = 0
intrachromosomal = True
notSkipped = 0 #just some counters to see how many regions have 0.0, these will be skipped
total = 0 #see how many lines are in the file in total
with open(subdir + regionsFile, 'w') as regionsOut:
	with open(subdir + relationsFile, 'w') as relationsOut:
		#Write the headers for the Neo4J input file
		regionsOut.write("id:ID\n")
		relationsOut.write(":START_ID,:END_ID\n")
		
			
		for hiCFilePath in filesToRead:
			hiCFile = glob.glob(hiCFilePath) #get the actual file path instead of the wildcard one
			print hiCFile[0] #There will always be only one file that matches
			
			#Read the chromosome identifier from the file name, this is not present in the data
			splitPath = re.split("/", hiCFile[0])
			fileName = splitPath[len(splitPath)-1]
		
			splitFileName = re.split("_", fileName)
		
			chr1Split = list(splitFileName[0])
			chr1 = "".join(chr1Split[3:len(chr1Split)])
			chr2 = splitFileName[1] 
		
			with open(hiCFile[0], "r") as f:
				for line in f:
					total += 1
					splitLine = line.split("\t")
					
					#Skip lines that have less than 1 interaction. 
					if splitLine[2] != "0.0": #This is faster than the larger than 0, when converting to float first
						notSkipped += 1
						start = "chr" + str(chr1) + "_" + splitLine[0]
						if intrachromosomal == True: #if intrachromosomal, chr2 is the same as chr1
							end = "chr" + str(chr1) + "_" + splitLine[1]
						else:	
							end = "chr" + str(chr2) + "_" + splitLine[1]
						
						#We wish to only write a region to the regions file (which will be a node in the graph), when it has not already been written to the file previously.
						#To save computational time, we add a key to the dictionary, dictionaries are really fast. The value does not matter.
						#If the length changes after adding the key, which is a very wuick check, we know that we have a new region. We write this to the file.
						#This is repeated for the end position of the interaction. 
						seenPositions[start] = 0 #First check if the start position changes the lenght
						
						diff = len(seenPositions) - previousLen #check if we added something or not. 
						if diff == 1:
							#write start to file
							regionsOut.write(start + "\n")
							previousLen = len(seenPositions)
							
						seenPositions[end] = 0
						 
						diff = len(seenPositions) - previousLen
						if diff == 1:
							#write end to file
							regionsOut.write(end + "\n")
							previousLen = len(seenPositions)
						
						relationsOut.write(start + "," + end + "\n")
					
				f.close()
				#break #if you wish to read only one file for testing, uncomment this

endTime = time.time()
print "Took ", endTime - startTime, " seconds to store edges from all files"

#Upload the data to Neo4J
os.system("./importData.sh")


