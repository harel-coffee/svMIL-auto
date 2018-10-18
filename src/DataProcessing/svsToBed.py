#As VCF does not work, convert to BED instead


import sys
import numpy as np

variantsFile = sys.argv[1]
outFile = sys.argv[2]

variantsList = []

with open(variantsFile, 'rb') as f:
	
	lineCount = 0
	header = []
	for line in f:
		line = line.strip()
		splitLine = line.split("\t")
		
		#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
		if lineCount < 1:

			header = splitLine
			lineCount += 1
			continue
		
		#Now extract the chromosome, start and end (there are multiple)
		chr1Index = header.index("chr1")
		s1Index = header.index("s1")
		e1Index = header.index("e1")

		chr2Index = header.index("chr2")
		s2Index = header.index("s2")
		e2Index = header.index("e2")

		cancerTypeIndex = header.index("cancer_type")
		sampleNameIndex = header.index("sample_name")
		svTypeIndex = header.index("sv_type")
		
		cancerType = splitLine[cancerTypeIndex]
		sampleName = splitLine[sampleNameIndex]
		svType = splitLine[svTypeIndex]
		#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
		if splitLine[chr1Index] == splitLine[chr2Index]:
			if splitLine[s2Index] == 'NaN':
				splitLine[s2Index] = int(splitLine[s1Index])
				
			if splitLine[e2Index] == 'NaN':
				splitLine[e2Index] = int(splitLine[e1Index])
		else:
			if splitLine[chr2Index] == 'NaN':
				continue # This line does not have correct chromosome 2 information (should we be skipping it?)

		s1 = str(splitLine[s1Index])
		e1 = str(splitLine[e1Index])
		s2 = str(splitLine[s2Index])
		e2 = str(splitLine[e2Index])
		chr2 = splitLine[chr2Index]
		
		chr1 = splitLine[chr1Index]
		
		#Some positions are swapped
		if int(e2) < int(e1):
			tmpE1 = e1
			e1 = e2
			e2 = tmpE1
			tmpS1 = s1
			s1 = s2
			s2 = tmpS1
		
		#Sometimes only the end is swapped.
		if int(e2) < int(s2):
			tmpS2 = s2
			s2 = e2
			e2 = tmpS2
			
		if int(e1) < int(s1):
			tmpS1 = s1
			s1 = e1
			e1 = tmpS1
			
			
		#chr 1, start, end, chr2, start2, end2
		variantsList.append([chr1, s1, e1, chr2, s2, e2, svType, cancerType, sampleName])

regions = np.array(variantsList, dtype='object')

#Write to bed format

with open(outFile, 'wb') as out:
	
	for region in regions:
		
		concatenatedPosition = region[0] + ":" + region[1] + "-" + str(int(region[2]) + 1) #Make sure to do the +1 because otherwise bedtools won't recognize it as one base
		
		line = "chr" + region[0] + "\t" + region[1] + "\t" + str(int(region[2]) + 1) + "\t" + concatenatedPosition + "\n"
		out.write(line)
		

#For translocations, we only need the 'main' chromosome. So for the end position we can always use the 'start' of the translocation, assuming chromosome 1.		
	



