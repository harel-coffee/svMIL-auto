
import sys
import numpy as np

burdenFile = sys.argv[1]
patternFile = sys.argv[2]
combineType = sys.argv[3] #multiply or correlation?
outFile = sys.argv[4]

burden = np.loadtxt(burdenFile, dtype="object")
pattern = np.loadtxt(patternFile, dtype="object")

combined = np.empty(burden.shape, dtype="object")
combined[:,0] = burden[:,0]
for row in range(0, burden.shape[0]):
	
	burdenScore = burden[row,burden.shape[1]-1]
	
	#find the right gene
	geneName = burden[row,0]
	
	genePattern = pattern[pattern[:,0] == geneName][0]
	patternScore = genePattern[pattern.shape[1]-1]
	
	if combineType == "Multiply":
		totalScore = float(burdenScore) * float(patternScore)
	
	if combineType == "Correlation":
		burdenScores = [float(burdenVal) for burdenVal in burden[row,1:27]]
		patternScores = [float(patternVal) for patternVal in pattern[row,1:27]]
		totalScore = np.corrcoef(burdenScores, patternScores)[0,1]
		
	
	combined[row, combined.shape[1]-1] = totalScore
	
burdenTotal = [float(burden) for burden in burden[:,28]]
patternTotal = [float(pattern) for pattern in pattern[:,28]]
print np.corrcoef(burdenTotal, patternTotal)


combined = combined[combined[:,28].argsort()[::-1]] #Select the column  to rank by	
		
header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\ttotal"
				
#Write to numpy output file	
np.savetxt(outFile, combined, delimiter='\t', fmt='%s', header=header)	


	

