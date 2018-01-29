#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""	
	Input: tab-delimited file with variants (or simply, regions). Required format: columns chr 1, s1, e1, chr2, s2, e2 per region on the rows. 
	Output: null
	Functionality: starting point of the annotation pipeline. Sorts the input file and starts the annotation process.
"""
### Imports ###
import sys
import numpy as np
from annotator import Annotator

### Code ###

def parseInputData():
	"""
		Parse the regional information from the input file provided on the command line. If s2 or e2 are NaN, the values for s1 and e1 are used instead and are copied to s2 and e2.
		Positions are converted to integers.
		
		Returns:
			inputData (list): list with regions, where each entry is a list with entries in the order: chr 1, s1, e1, chr2, s2, e2
	"""
	
	inFile = sys.argv[1]
	
	inputData = []
	with open(inFile, "r") as f:
		lineCount = 0
		for line in f:
			if lineCount < 2:
				lineCount += 1
				continue
			line = line.strip()
			splitLine = line.split("\t")
			
			#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
			if splitLine[0] == splitLine[4]:
				if splitLine[5] == 'NaN':
					splitLine[5] = int(splitLine[1])
					
				if splitLine[6] == 'NaN':
					splitLine[6] = int(splitLine[2])
			else:
				if splitLine[5] == 'NaN':
					continue # This line does not have correct chromosome 2 information
			
			#Convert the start and end positions to int here so that we do not need to do this in the loop, it slows down the code
			splitLine[1] = int(splitLine[1])
			splitLine[2] = int(splitLine[2])
			splitLine[5] = int(splitLine[5])
			splitLine[6] = int(splitLine[6])
			
			#chr 1, start, end, chr2, start2, end2
			inputData.append([splitLine[0], splitLine[1], splitLine[2], splitLine[4], splitLine[5], splitLine[6]])
	
	
	
	return inputData

def prepareInputData(inputData):
	"""
		Converts the input data to a Numpy matrix for speed, and then sorts it (chr1 < chr2, s1 < e1, s2 < e2).
		
		Parameters:
			inputData (list): output of function parseInputData
			
		Returns:
			sortedInputData (numpy matrix): inputData sorted (chr1 < chr2, s1 < e1, s2 < e2)
	"""
	#Convert the pancancerData to numpy format so that we can do stuff like sorting later
	inputDataNp = np.array(inputData)
	
	chr1col = inputDataNp[:,0]
	chr2col = inputDataNp[:,3]
	
	sortedInd = np.lexsort((chr2col, chr1col)) #sort first by column 1, then by column 2. This works, but it is lexographical, so chromosome 11 comes before chromosome 2. For this purpose it is ok, since we only look for this
												#chromosome in other files to subset these, so the order in which we do that does not really matter. 
	
	sortedInputData = inputDataNp[sortedInd]
	
	return sortedInputData


#1. Initialize
annotator = Annotator()

#2. Read input data
inputData = parseInputData()

#3. Prepare input data, including sorting for speed
preparedInputData = prepareInputData(inputData)

#4. Do the annotation
annotator.annotate(preparedInputData)

