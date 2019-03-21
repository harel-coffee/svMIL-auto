import sys
import numpy as np

#1. Get the total number of samples and choose 89 of these at random

uniqueSamples = []
header = dict()
with open(sys.argv[1], 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			
			for entryInd in range(0, len(splitLine)):
				header[splitLine[entryInd]] = entryInd
			lineCount += 1
			continue
		
		if len(splitLine) == header['samples']: #if there are no samples associated to this variant (which is weird tbh???)
			continue
	
		samples = splitLine[header['samples']]
		
		splitSamples = samples.split(",")
		
		for sample in splitSamples:
			if sample not in uniqueSamples:
				uniqueSamples.append(sample)

print uniqueSamples		
print len(uniqueSamples)

sampleSubset = np.random.choice(uniqueSamples, 89, replace=False)

#2. Collect all SVs from these samples, check if the number is about equal to the somatic SVs
subsetSVs = []
with open(sys.argv[1], 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			lineCount += 1
			continue
		
		if len(splitLine) == header['samples']: #if there are no samples associated to this variant (which is weird tbh???)
			continue
		
		if splitLine[5] == "deletion" or splitLine[5] == "duplication" or splitLine[5] == "inversion":

			samples = splitLine[header['samples']]
			
			splitSamples = samples.split(",")
			
			for sample in splitSamples:
				if sample in sampleSubset:
					
					chr1 = splitLine[1]
					s1 = splitLine[2]
					e1 = splitLine[3]
					
					o1 = "+"
					o2 = "-"
					
					source = "dgv"
					sample_name = sample
					sv_type = splitLine[5]
					cancer_type = "germline"
					
					subsetSVs.append([chr1, s1, e1, o1, chr1, s1, e1, o2, source, sample_name, sv_type, cancer_type])

print len(subsetSVs)

#Subsample within the SV subset to get a similar number of SVs as for the somatic case

somaticSvNum = 17649
subsetLines = np.random.choice(range(0, len(subsetSVs)), somaticSvNum, replace=False)

subsampledSVs = []

lineCount = 0
with open(sys.argv[2], 'w') as outF:
	outF.write("chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type\n")
	for sv in subsetSVs:
		
		if lineCount in subsetLines:
			outF.write("\t".join(sv))
			outF.write("\n")

			subsampledSVs.append(sv)
			
		lineCount += 1
			
print len(subsampledSVs)	




