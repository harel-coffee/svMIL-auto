import sys
import re
import numpy as np

filteredSvs = []

with open(sys.argv[1], 'r') as inF:
		
		for line in inF:
			if re.match("^#", line):
				continue
			
			splitLine = line.split("\t")
			
			if splitLine[6] == "PASS":
				if splitLine[4] == "DUP" or splitLine[4] == "DEL" or splitLine[4] == "INV":
					
					#Output to the right format chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type
					chr1 = splitLine[0]
					s1 = splitLine[1]
					e1 = splitLine[2]
					
					o1 = "+"
					o2 = "-"
					
					source = "gnomad"
					sample_name = splitLine[3]
					sv_type = splitLine[4]
					cancer_type = "Germline"
					
					filteredSvs.append([chr1, s1, e1, o1, chr1, s1, e1, o2, source, sample_name, sv_type, cancer_type])

filteredSvs = np.array(filteredSvs, dtype="object")

somaticSvNum = 17649
#First determine how many SVs there are in total
#Then randomly sample numbers from this range
#When writing SVs to file, then we only write these random lines
svNum = filteredSvs.shape[0]

randomSvNums = np.random.choice(range(1,svNum), somaticSvNum, replace=False)

lineCount = 0
with open(sys.argv[2], 'w') as outF:
	outF.write("chr1	s1	e1	o1	chr2	s2	e2	o2	source	sample_name	sv_type	cancer_type\n")
	for sv in filteredSvs:
		
		if lineCount not in randomSvNums:
			lineCount += 1
			continue
		
		svLine = "\t".join(sv)
		outF.write(svLine)
		outF.write("\n")
		lineCount += 1