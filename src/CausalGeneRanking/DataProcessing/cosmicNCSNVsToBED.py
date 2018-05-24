
"""
	The purpose of this script is to convert the cosmic NC SNV file to a BED file that can be loaded into IGV.
	
	Format:
	
	chromosome	start	end	sampleName

"""

import sys


inFile = sys.argv[1]
outFile = sys.argv[2]


with open(outFile, 'wb') as outF:

	with open(inFile, 'rb') as f:
		lineCount = 0
		header = []
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")	
		
			if lineCount < 1: #get the header
				header = splitLine
				lineCount += 1
				continue
			
			positionInd = header.index("genome position") #chromosome and position
			tumorIDInd = header.index("Primary site")
			sampleIDInd = header.index("Sample name")
			
			position = splitLine[positionInd]
			tumorID = splitLine[tumorIDInd]
			
			sampleID = splitLine[sampleIDInd]
			
			
			if tumorID != "breast":
				continue #Make this file only for breast cancer for now. 
			
			print sampleID
			
			colonSplitPosition = position.split(":")
			dashSplitPosition = colonSplitPosition[1].split("-")
			
			chromosome = colonSplitPosition[0]
			start = dashSplitPosition[0]
			end = dashSplitPosition[1]
			
			newLine = chromosome + "\t" + start + "\t" + end + "\t" + sampleID + "\n"
			outF.write(newLine)
				
			

	


