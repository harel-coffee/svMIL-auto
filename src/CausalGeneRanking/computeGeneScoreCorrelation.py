"""
	Output a file where the total score is the correlation between the scores of the individual features for each gene. 

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from six.moves import range

scoreFile = sys.argv[1]
outfileName = sys.argv[2]

scores = np.loadtxt(scoreFile, dtype="object")
print(scores)

correlation = np.empty(scores.shape, dtype="object")

for row in range(0, scores.shape[0]):
	
	correlation[row,0:scores.shape[1]-1] = scores[row,0:scores.shape[1]-1] 
	
	
	#Get the gene scores but without the total and without the gene name
	geneScores = [float(score) for score in scores[row,1:scores.shape[1]-1]]
	
	#Correlate these scores
	corr = np.corrcoef(geneScores)
	correlation[row,scores.shape[1]-1] = corr
	
	
#Write the correlation to a new file
header = "geneName\tgeneScore\teQTLGains\teQTLLosses\tenhancerGains\tenhancerLosses\tpromoterGains\tpromoterLosses\tcpgGains\tcpgLosses\ttfGains\ttfLosses\thicGains\thicLosses\th3k9me3Gains\th3k9me3Losses\th3k4me3Gains\th3k4me3Losses\th3k27acGains\th3k27acLosses\th3k27me3Gains\th3k27me3Losses\th3k4me1Gains\th3k4me1Losses\th3k36me3Gains\th3k36me3Losses\tdnaseIGains\tdnaseILosses\tcorrelation"
				
#Write to numpy output file	
np.savetxt(outfileName, correlation, delimiter='\t', fmt='%s', header=header)

	
	
	
