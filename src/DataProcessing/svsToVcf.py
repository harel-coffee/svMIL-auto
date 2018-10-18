#Read SVs of a specific tumor type, then write these to a VCF file

#Or just write to a BED file and see how it works with multiple samples.  
import sys
import numpy as np

variantsFile = sys.argv[1]
outFile = sys.argv[3]
refFile = sys.argv[2] #breastGastricSvs.vcf

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
		
		#Sometimes chr1 and chr2 are swapped
		if int(e2) < int(e1):
			tmpE1 = e1
			e1 = e2
			e2 = tmpE1
			tmpS1 = s1
			s1 = s2
			s2 = tmpS1
	
		#chr 1, start, end, chr2, start2, end2
		variantsList.append([chr1, s1, e1, chr2, s2, e2, svType, cancerType, sampleName])

regions = np.array(variantsList, dtype='object')

#Filter by cancer type
cancerType = "breast/gastric"
filteredRegions = regions[np.where(regions[:,7] == cancerType),:][0]

#Read the REF information by coordinate
#This is a vcf file

refBases = dict()

with open(refFile, 'r') as ref:
	
	lineCount = 0
	
	for line in ref:
		# 
		# if lineCount < 2:
		# 	lineCount += 1
		# 	continue
		# 
		# line = line.strip()
		# 
		# splitLine = line.split("\t")
		# 
		# #coordinates key is made up out of the chrom:start
		# chrom = splitLine[0]
		# 
		# info = splitLine[7]
		# splitInfo = info.split(";")
		# 
		# e1 = splitInfo[3]
		# splitE1 = e1.split("=")
		# endPos = splitE1[1]
		# 
		# 
		# #Here also, if the positions are switched, use the end as pos instead of start.
		# if int(endPos) < int(splitLine[1]):
		# 	pos = endPos
		# else:
		# 	pos = splitLine[1]
		# 
		# coordinates = chrom + ":" + pos
		# 
		# refBase = splitLine[3]
		# 
		# refBases[coordinates] = refBase
		
		coordinates = splitLine[0]
		refBase = splitLine[1]
		refBases[coordinates] = refBase


#Write all these to a VCF file, add the sample name

with open(outFile, 'wb') as out:
	
	out.write("##fileformat=VCFv4.0\n")
	out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
	
	for region in filteredRegions:
		#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
		#INFO field should contain the SVTYPE, and also the chr2 information for translocations
		
		#The REF information should be parsed from the file generated with bedtools containing the reference bases at that position. 

		coordinates = region[0] + ":" + region[1] + "-" + str(int(region[2]))
		#coordinates = region[0] + ":" + region[1]
		
		end = region[5]
		if region[0] != region[3]:
			end = region[2]
		
		refBase = refBases[coordinates][0] #only show the first base, otherwise there is way too much data to show in IGV, some SVs are very long
		#refBase = refBases[coordinates]
		
		sampleName = region[8].replace(" ", "")
		info = "SVTYPE=" + region[6] + ";CHR2=" + region[3] + ";S2=" + region[4] + ";E1=" + region[2] + ";SAMPLE=" + sampleName + ";END=" + end #IGV uses 'END' for plotting, so don't use the e2 coordinate if translocation
		
		line = region[0] + "\t" + region[1] + "\t.\t" + refBase + "\t.\t.\t.\t" + info + "\n"
		
		out.write(line)
		
		
		
