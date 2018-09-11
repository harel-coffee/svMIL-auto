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
		- Clarity of how paths are being read

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

#hiCFolder = "../../../../../Data/HiC/HUVEC_interchromosomal/5kb_resolution_interchromosomal/" #Use and replace this path if you want to run the script yourself with your own local dataset
hiCFolder = "../../../data/thresholdedHiC/" #Edit to read the thresholded data. The code will read a folder inside this thresholdedHiC folder, so if the folder Intrachromosomal is in there, it will read from inside that Intrachromosomal directory. 

#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
filesToRead = []
for folderName in os.listdir(hiCFolder):
	
	#Skip if it is an osx specific file and not a folder
	if os.path.isdir(hiCFolder + "/" + folderName):
		filesToRead = glob.glob(hiCFolder + "/" + folderName + "/*.RAWobserved_thresholded_5")

#Read the hg19 genes
hg19Genes = "../../../data/Genes/ensemblGenesHg19"

geneLocations = []
with open(hg19Genes, 'r') as geneIn:
	lineCount = 0
	for line in geneIn:
		line = line.strip()
		
		if lineCount < 1:
			lineCount += 1
			continue
		splitLine = line.split("\t")
		
		#chr, start, end, name
		geneLocations.append([splitLine[1], int(splitLine[3]), int(splitLine[4]), splitLine[6]])

geneLocations = np.array(geneLocations, dtype='object')

#Read the causal genes
causalGenes = "../../../data/Genes/causalGenes.bed"
causalGeneNames = dict() #easy searching if a gene is causal yes/no
causalGeneLocations = []
with open(causalGenes, 'r') as geneIn:
	lineCount = 0
	for line in geneIn:
		line = line.strip()
		
		if lineCount < 1:
			lineCount += 1
			continue
		splitLine = line.split("\t")
		
		if splitLine[3] not in causalGeneNames:
			causalGeneNames[splitLine[3]] = 0
			
		#chr, start, end, name
		causalGeneLocations.append([('chr' + splitLine[0]), int(splitLine[1]), int(splitLine[2]), splitLine[3]])

causalGeneLocations = np.array(causalGeneLocations, dtype='object')
	
startTime = time.time()

#The output files to write the results to, these will be input for Neo4J
regionsFile = "regions.csv"
relationsFile = "regions_regions_rel.csv"

#This subdir is where Neo4J will write the data to
subdir = "Intrachromosomal/"

if not os.path.exists(subdir):
   os.makedirs(subdir)


seenPositions = dict() #Keep a dictionary with previous regions that were already seen. Neo4J has issues with resolving collisions when there are many non-unique regions, so filtering here saves computational time. 

previousLen = 0
intrachromosomal = True #Used to see if chr1 and chr2 need to be parsed, or only chr1
notSkipped = 0 #just some counters to see how many regions have 0.0, these will be skipped
total = 0 #see how many lines are in the file in total
with open(subdir + regionsFile, 'w') as regionsOut:
	with open(subdir + relationsFile, 'w') as relationsOut:
		#Write the headers for the Neo4J input file
		regionsOut.write("id:ID,geneNames:string,causality:string\n")
		relationsOut.write(":START_ID,:END_ID\n")
		
		for hiCFile in filesToRead:
			print hiCFile
			#Read the chromosome identifier from the file name, this is not present in the data
			splitPath = re.split("/", hiCFile)
			fileName = splitPath[len(splitPath)-1]
		
			splitFileName = re.split("_", fileName)
		
			chr1Split = list(splitFileName[0])
			chr1 = "".join(chr1Split[3:len(chr1Split)])
			chr2 = splitFileName[1] 
		
			with open(hiCFile, "r") as f:
				for line in f:
					total += 1
					splitLine = line.split("\t")
					
					#Skip lines that have less than 1 interaction 
					if splitLine[2] != "0.0": #This is faster than the larger than 0, when converting to float first
						notSkipped += 1
						start = "chr" + str(chr1) + "_" + splitLine[0]
						if intrachromosomal == True: #if intrachromosomal, chr2 is the same as chr1
							end = "chr" + str(chr1) + "_" + splitLine[1]
						else:	
							end = "chr" + str(chr2) + "_" + splitLine[1]
						if start == end:
							continue #avoid adding self-loops to the database, these are not useful for now
							
						
						
						#We wish to only write a region to the regions file (which will be a node in the graph), when it has not already been written to the file previously.
						#To save computational time, we add a key to the dictionary, dictionaries are really fast. The value does not matter.
						#If the length changes after adding the key, which is a very wuick check, we know that we have a new region. We write this to the file.
						#This is repeated for the end position of the interaction. 
						seenPositions[start] = 0 #First check if the start position changes the lenght
						
						diff = len(seenPositions) - previousLen #check if we added something or not. 
						if diff == 1:
							#write start to file
							regionsOut.write(start + ",")
							
							#Check if this region overlaps with a gene or causal gene
							#First get the subset with the right chromosome
							chrSubset = geneLocations[geneLocations[:,0] == 'chr' + str(chr1)]
							
							
							startEnd = int(splitLine[0]) + 5000 #the end coordinate of the region starting the interaction
							
							regionStartMatches = int(splitLine[0]) <= chrSubset[:,2]
							regionEndMatches = startEnd >= chrSubset[:,1]

							allMatches = regionStartMatches * regionEndMatches
							
							matchingGenes = chrSubset[allMatches,:] #the genes in this region, can be multiple
							
							geneNames = dict()
							causality = False
							if len(matchingGenes) > 0:
								tmpGeneNames = matchingGenes[:,3]
								for geneName in tmpGeneNames: #use a dictionary to quickly make genes unique
									
									if geneName in causalGeneNames:
										causality = True
									geneNames[geneName] = 0
								
							geneNames = ";".join(geneNames)
							
								
							regionsOut.write(geneNames + ",")
							regionsOut.write(str(causality) + "\n")
							previousLen = len(seenPositions)
							
						seenPositions[end] = 0
						 
						diff = len(seenPositions) - previousLen
						if diff == 1:
							#write end to file
							regionsOut.write(end + ",")
							
							
							#Check if this region overlaps with a gene or causal gene
							#First get the subset with the right chromosome
							if intrachromosomal == True:
								chrSubset = geneLocations[geneLocations[:,0] == 'chr' + str(chr1)]
							else:
								chrSubset = geneLocations[geneLocations[:,0] == 'chr' + str(chr2)]
							
							startEnd = int(splitLine[1]) + 5000 #the end coordinate of the region ending the interaction
							
							regionStartMatches = int(splitLine[1]) <= chrSubset[:,2]
							regionEndMatches = startEnd >= chrSubset[:,1]
							
							allMatches = regionStartMatches * regionEndMatches
							
							matchingGenes = chrSubset[allMatches] #the genes in this region, can be multiple
							
							geneNames = dict()
							causality = False
							if len(matchingGenes) > 0:
								tmpGeneNames = matchingGenes[:,3]
								for geneName in tmpGeneNames: #use a dictionary to quickly make genes unique
									
									if geneName in causalGeneNames:
										causality = True
								
									geneNames[geneName] = 0
								
							geneNames = ";".join(geneNames)
							
							regionsOut.write(geneNames + ",")
							regionsOut.write(str(causality) + "\n")

							previousLen = len(seenPositions)
						
						relationsOut.write(start + "," + end + "\n")
					
				f.close()
				#break #if you wish to read only one file for testing, uncomment this

endTime = time.time()
print "Took ", endTime - startTime, " seconds to store edges from all files"

#Upload the data to Neo4J
if intrachromosomal == True:
	os.system("./importData_intrachromosomal.sh")
else:
	os.system("./importData_interchromosomal.sh")

