"""
	Purpose: read a .seg file from the TCGA data, and convert the segment mean values to copy numbers.
	From these data, we can determine where deletions are located in the patients. 

"""

import sys
from math import pow
import numpy as np

inFile = sys.argv[1]
delOutFile = sys.argv[2]

copyNumbers = []

with open(inFile, 'r') as inF:
	lineCount = 0
	for line in inF:
		
		if lineCount < 1:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		segMean = float(splitLine[5])
		
		
		#(2^seg_mean)*2
		
		cn = pow(2,segMean)*2
		
		chrom = splitLine[1]
		if splitLine[1] == "23":
			chrom = "X"
		
		#sampleID, chrom, start, end, segMean, cn
		copyNumbers.append([splitLine[0], chrom, splitLine[2], splitLine[3], segMean, cn])

copyNumbers = np.array(copyNumbers, dtype='object')

print copyNumbers

#write only the deletions to a file

with open(delOutFile, 'w') as outF:
	outF.write("chr1\ts1\te1\tchr2\ts2\te2\tcancer_type\tsample_name\tsv_type\n") #write a header
	
	for cn in copyNumbers:
		
		
		if cn[5] < 2: #everything below the standard of 2 can be considered a deletion (in at least one subclone of the tumor)
			
			outF.write(cn[1] + "\t" + cn[2] + "\t" + cn[3] + "\t" + cn[1] + "\t" + cn[2] + "\t" + cn[3] + "\t" + "breast\t" + cn[0] + "\tdeletion\n")
			
			
		
	


