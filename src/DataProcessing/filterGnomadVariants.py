"""
	This script parses SVs from gnomAD as an alternative to DGV.

	Because there are only 8 translocations in gnomAD, there is not enough power to
	say something about what those translocations affect in germline more than by
	random chance. Therefore we only extract DEL/DUP/INV SVs.

	SVs are subsampled to the same number as the number of SVs in the HMF BRCA dataset.
	SV types are not matched, as DELs are enriched for germline, which may hold
	important information about the distribution & TAD effects in germline.

"""


import sys
import numpy as np
np.random.seed(785)


#Get the variants from gnomAD and parse these to a format that we can use
subsetSVs = []
with open(sys.argv[1], 'r') as inF:
	lineCount = 0
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		if lineCount < 1:
			lineCount += 1
			continue

		if splitLine[4] == "DEL" or splitLine[4] == "DUP" or splitLine[4] == "INV":

			sample = splitLine[3] #there does not seem to be sample information so use the ID to keep it unique
			#for germline it does not matter since we only use it for the heatmap, where grouping by sample is not necessary.

			chr1 = 'chr' + splitLine[0]
			s1 = splitLine[1]
			e1 = splitLine[2]

			#dummy value
			o1 = "+"
			o2 = "-"

			source = "gnomAD"
			sample_name = sample
			sv_type = splitLine[4]
			cancer_type = "germline"

			subsetSVs.append([chr1, s1, e1, o1, chr1, s1, e1, o2, source, sample_name, sv_type, cancer_type])

#Subsample within the SV subset to get a similar number of SVs as for the somatic case

somaticSvNum = 73293
np.random.seed(785)
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




