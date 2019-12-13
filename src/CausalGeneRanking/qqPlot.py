
import numpy as np
import matplotlib.pyplot as plt

#Load the z-score ranks for both cases
tadRanks = np.loadtxt('pValues_ranks.txt', dtype='object')
ruleRanks = 0


plt.scatter(np.sort(tadRanks[:,6]), np.sort(ruleRanks))
plt.xlabel('Ranks for TAD-based')
plt.ylabel('Ranks for rule-based')
plt.show()