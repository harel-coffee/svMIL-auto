import numpy as np


ranks = np.loadtxt('tadBasedRanks.rnk', dtype='object')

ranks[:,1] = ranks[:,1].astype(float)

print(np.max(ranks[:,1]))
print(np.min(ranks[:,1]))
	
	
for rank in ranks:
	
	if float(rank[1]) > 180:
		print(rank)
