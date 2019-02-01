#Take ensembl gff gene file as input
#Return file with only the genes

import re
import sys
import numpy as np

inFile = sys.argv[1]
outFile = sys.argv[2]

genes = dict()
with open(inFile) as inF:
	
	for line in inF:
		
		if re.search("#", line) == None:
			line = line.strip()
			splitLine = line.split("\t")
			
			if splitLine[2] == "gene":
				
				info = splitLine[8]
				splitInfo = info.split(";")
				fieldDict = dict()
				for field in splitInfo:
					splitField = field.split("=")
					
					fieldDict[splitField[0]] = splitField[1]
				
				if fieldDict["gene_biotype"] != "protein_coding":
					continue
				# 	
				# if splitField[0] == "Name":
				# 	
				# 	if splitField[1] == "LOC105377838":
				# 		print "LOC105377838"
				# 		print line
				# 		exit()
					
				if fieldDict["Name"] not in genes:
					genes[fieldDict["Name"]] = []
				#Complicated chromosome encoding
				splitLocusId = splitLine[0].split(".")
				
				splitSplitLocusId = splitLocusId[0].split("_")
				
				#Now replace all 0s in the second part
				#But make sure that numbers such as 10 or 20 are not removed
				#chromId = splitSplitLocusId[1].replace("0", "")
				# print splitSplitLocusId[1]
				chromId = list(splitSplitLocusId[1])
				
				numberFound = False
				chromNumber = ""
				for char in chromId:
					if char != "0":
						numberFound = True
					if numberFound == True:	
						chromNumber += char	
					
				
				# print chromId
				# print chromNumber
				# exit()
				
				
				
				genes[fieldDict["Name"]].append([fieldDict["Name"], 'chr' + chromNumber, splitLine[3], splitLine[4]])
				
			
		
			
				
				
			
				
filteredGenes = []			
for gene in genes:
	
	if len(genes[gene]) > 1: #if there are more genes with the same name, get the longest
		currentLongest = genes[gene][0]
		currentLength = int(genes[gene][0][3]) - int(genes[gene][0][2])
		for transcript in genes[gene]:
			length = int(transcript[3]) - int(transcript[2])
			if length > currentLength:
				currentLength = length
				currentLongest = transcript
		
		filteredGenes.append(currentLongest)
	else:
		filteredGenes.append(genes[gene][0])
		
print len(filteredGenes)
print np.array(filteredGenes)

#Write the genes to an output file

with open(outFile, 'w') as outF:
	
	for gene in filteredGenes:
		outF.write(gene[1] + "\t" + gene[2] + "\t" + gene[3] + "\t" + gene[0] + "\n")
		
	
		


