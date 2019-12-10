"""
	Get the DEG set from the TAD disruption analysis and the one from linking SVs to genes. Which overlap? Are there things missing?

"""

import sys
import numpy as np

tadDegs = np.loadtxt('tadDisr/pValues_shuffled_0.txt', dtype='object')
svDegs = np.load(sys.argv[1], allow_pickle=True, encoding='latin1')

degs = 0
for deg in tadDegs:
	
	if deg[3] == 'True':
		degs += 1
		print(deg)

print(degs)
exit()

print(svDegs.shape)
print(tadDegs)

missingSets = []
for deg in svDegs:
	
	splitDeg = deg[0].split('_')
	
	gene = splitDeg[0]
	patient = splitDeg[7]
	
	combi = patient + '_' + gene
	
	
	
	if combi in tadDegs:
		
		print('overlap: ', combi)
		match = tadDegs[tadDegs[:,0] == combi]
		print(match)
		missingSets.append(combi)
	
#np.savetxt('missingSets.txt', missingSets, fmt='%s')	
	

