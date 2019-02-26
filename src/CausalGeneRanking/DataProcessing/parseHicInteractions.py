"""
	Goal: read the thresholded Hi-C interactions and write these to a singular file with interactions
	
	1. Filter for interactions that start/end within one TAD.
	2. Filter for interactions that take place with a bin containing a gene
	
"""

import os
import glob
import re

hiCFolder = "../../../data/thresholdedHiC/" #Edit to read the thresholded data. The code will read a folder inside this thresholdedHiC folder, so if the folder Intrachromosomal is in there, it will read from inside that Intrachromosomal directory. 
intrachromosomal = True
#Obtain the names of the files that we wish to read for each folder, use a wildcard since the names are different
filesToRead = []
for folderName in os.listdir(hiCFolder):
	
	#Skip if it is an osx specific file and not a folder
	if os.path.isdir(hiCFolder + "/" + folderName):
		filesToRead = glob.glob(hiCFolder + "/" + folderName + "/*.RAWobserved_thresholded_5")

#1. Get the TADs


relationsFile = '../../../data/HiC/HMEC_interactions.txt'
with open(relationsFile, 'w') as relationsOut:
	
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
				splitLine = line.split("\t")
				
				#Skip lines that have less than 1 interaction 
				if splitLine[2] != "0.0": #This is faster than the larger than 0, when converting to float first
					start = "chr" + str(chr1) + "_" + splitLine[0]
					if intrachromosomal == True: #if intrachromosomal, chr2 is the same as chr1
						end = "chr" + str(chr1) + "_" + splitLine[1]
					else:	
						end = "chr" + str(chr2) + "_" + splitLine[1]
					if start == end:
						continue #avoid adding self-loops to the database, these are not useful for now
						

					relationsOut.write(start + "\t" + end + "\n")
