#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""	
	Input: tab-delimited file with variants (or simply, regions). Required format: columns chr 1, s1, e1, chr2, s2, e2, somatic (yes/no) per region on the rows. 
	Output: null
	Functionality: starting point of the annotation pipeline. Sorts the input file and starts the annotation process.
	
	TODO:
		- Update input data format
"""

### Imports ###
import sys
import numpy as np
from annotator import Annotator

### Code ###

def parseInputData():
	"""
		Parse the regional information from the input file provided on the command line. 
		Positions are converted to integers.
		
		For SVs:
			- If s2 or e2 are NaN, the values for s1 and e1 are used instead and are copied to s2 and e2.
		For SNVs/SNPs:
			- SNV/SNP files can be provided with only chr1, s1 and e1. For SVs, the values of s1 and e2 are used when chr1 and chr2 are the same. So in the case of SNVs/SNPs, we can copy the values of chr1 to chr2,
			  e1 to e2, s1 to e1 and e2 to s2.
		
		Returns:
			inputData (list): list with regions, where each entry is a list with entries in the order: chr 1, s1, e1, chr2, s2, e2, regionType. A a regionType is passed to the rest of the code as part of every
			line in the inputData list, which can be SV, SNV or SNP. 
			
	"""
	
	inFile = sys.argv[1]
	
	inputData = []
	with open(inFile, "r") as f:
		lineCount = 0
		header = []
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			if lineCount < 1:
				header = splitLine
				lineCount += 1
				continue
			
			#First, we check if the header contains SV-related values or not. If not, it is an SNV/SNP.
			regionType = "SNV"
			if 'chr2' in header: 
				regionType = "SV"
			else:
				somaticColumnIndex = header.index("somatic") #find which column describes if the variant is somatic or not
				#Find the value of the 'somatic' column
				if splitLine[somaticColumnIndex] == "no":
					regionType = "SNP"
			
			#Obtain required data for each data type. Search which column has this information
			chr1Index = header.index("chr1")
			s1Index = header.index("s1")
			e1Index = header.index("e1")
			
			if regionType == "SV": #format the data for SVs specifically. Skip the orientation information for now. We may need it later

				chr2Index = header.index("chr2")
				s2Index = header.index("s2")
				e2Index = header.index("e2")

				#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
				if splitLine[chr1Index] == splitLine[chr2Index]:
					if splitLine[s2Index] == 'NaN':
						splitLine[s2Index] = int(splitLine[s1Index])
						
					if splitLine[e2Index] == 'NaN':
						splitLine[e2Index] = int(splitLine[e1Index])
				else:
					if splitLine[chr2Index] == 'NaN':
						continue # This line does not have correct chromosome 2 information (should we be skipping it?)
			
				s1 = int(splitLine[s1Index])
				e1 = int(splitLine[e1Index])
				s2 = int(splitLine[s2Index])
				e2 = int(splitLine[e2Index])
				chr2 = int(splitLine[chr2Index])
			
			else: #make sure that the end position is in e2 and the start in s1. Copy chr1 to chr2
				
				s1 = int(splitLine[s1Index])
				e1 = int(splitLine[s1Index])
				s2 = int(splitLine[e1Index])
				e2 = int(splitLine[e1Index])
				chr2 = int(splitLine[chr1Index])
			
				
			chr1 = int(splitLine[chr1Index])
			
			#chr 1, start, end, chr2, start2, end2, regionType
			inputData.append([chr1, s1, e1, chr2, s2, e2, regionType])
	
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

