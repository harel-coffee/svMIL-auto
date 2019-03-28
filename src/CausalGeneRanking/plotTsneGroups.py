"""
	Plot the tSNE for 2 groups of SVs.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt

#First, use non-coding scores > 0 vs = 0. Any specific SVs showing up that are different from the rest?

geneScores = np.loadtxt(sys.argv[1], dtype="object")

totalScores = np.array(geneScores[:,4:30], dtype="float")


# from sklearn.decomposition import PCA
# 
# pca = PCA(n_components=2)
# 
# projected = pca.fit_transform(totalScores)
# 
# colorLabels = []
# for score in geneScores:
# 	if float(score[30]) > 0:
# 		colorLabels.append('r')
# 	else:
# 		colorLabels.append('b')
# 
# plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
# plt.show()
# exit()
#Do a tSNE

from tsne import bh_sne


colorLabels = []
for score in geneScores:
	if float(score[30]) > 0:
		colorLabels.append('r')
	else:
		colorLabels.append('b')

vis_data = bh_sne(totalScores)
print vis_data
exit()
# plot the result
vis_x = vis_data[:, 0]
vis_y = vis_data[:, 1]



plt.scatter(vis_x, vis_y, c=colorLabels)
plt.show()
exit()
# 
# print labelList
# exit()
# 
# plt.scatter(vis_x, vis_y, c=labelList, edgecolor = 'none', alpha = 0.5, cmap=plt.cm.get_cmap("jet", 2))
# plt.colorbar()
# plt.show()
# 
# exit()
# 
# 
