#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

"""
	This script reads the original TP.txt and TN.txt files for SVs and converts them to an input format that can be recognized by the annotation pipeline prototype.
	All we need to do is update the header and make sure that there are no additional lines of text or ## signs in the header. 
	
	The accepted input format is (containing AT LEAST these fields):
	chr1	s1	e1	chr2	s2	e2	somatic
	
	The somatic column will have a yes/no value to distinguish SNPs from SNVs.
	The order of the columns is not important, and there can be multiple additional columns, which will simply not be used yet. 

"""

### Imports ###

import sys

### Code ###

#1. Read the file
inFile = sys.argv[1]
outFile = sys.argv[2]

with open(outFile, "w") as outF:
	
	with open(inFile, "r") as f:
		lineCount = 0
		header = ""
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			if lineCount < 1: #skip the text line in the header
				lineCount += 1
				continue
			#2. Update the header, make sure that the file starts with a header and not with text. Also remove the "##" from the start of the header, this is not necessary. 
			if lineCount < 2:
				splitHeader = splitLine
				#Remove the '##' from the chr1 value
				removedNumSign = [splitHeader[0][i] for i in range(2, len(splitHeader[0]))]
				removedNumSign = "".join(removedNumSign) #stitch the array back together into a string
				splitHeader[0] = removedNumSign #replace with the real chr1
				header = "\t".join(splitHeader) #join the header back together.
				
				#Write the header to the new file
				outF.write(header)
				outF.write("\n")
				
				lineCount += 1
				continue
			
			outF.write(line)
			outF.write("\n")
			lineCount += 1





