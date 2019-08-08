"""
	Make plots of the AUC and accuracy of the classifiers. 

"""

from __future__ import absolute_import
import matplotlib.pyplot as plt
import numpy as np
from six.moves import range

#For > 2 DEG labels
alphas = ['1e-15', '1e-10', '1e-8', '1e-4', '1e-3', '1e-2', '1', '5', '10', '20']
meanAcc = [0.15342345847111455, 0.15342345847111455, 0.15342345847111455, 0.15342345847111455, 0.15359498505773547, 0.15565448297407056, 0.15616965217237536, 0.15616965217237536, 0.15616965217237536, 0.15616965217237536]
meanAuc = [0.449106674438, 0.449106674438, 0.449106674438, 0.449018504241, 0.448981916369, 0.450752989686, 0.578084826086, 0.578084826086, 0.578084826086, 0.578084826086]
conceptPairs = [81996.9, 81921.2, 79934.3, 1544.8, 733.9, 1023.2, 0, 0,0,0]

fig, ax = plt.subplots(1)
plt.bar(np.arange(len(alphas)), meanAcc)
index = list(range(len(alphas)))
plt.xticks(index, alphas)
plt.show()

exit()

# fig, ax = plt.subplots(1)
# plt.bar(np.arange(len(alphas)), meanAuc)
# index = range(len(alphas))
# plt.xticks(index, alphas)
# #ax.set_xticklabels(alphas)
# plt.show()
# 
# fig, ax = plt.subplots(1)
# plt.bar(np.arange(len(alphas)), conceptPairs)
# index = range(len(alphas))
# plt.xticks(index, alphas)
# plt.show()


#Per-patient labels
alphas = ['1e-15', '1e-10', '1e-8', '1e-4', '1e-3', '1e-2', '1', '5', '10', '20']
meanAcc = [0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975, 0.933414047847975]
meanAuc = [0.969781504759, 0.969781504759, 0.969781504759, 0.970258732997, 0.970192873678, 0.968823113412, 0.966707023924, 0.966707023924, 0.966707023924, 0.966707023924]
conceptPairs = [87276.1, 87264.9, 84471.0, 859.4, 668.2, 179.5, 0, 0, 0, 0]

fig, ax = plt.subplots(1)
plt.bar(np.arange(len(alphas)), meanAcc)
index = list(range(len(alphas)))
plt.xticks(index, alphas)
plt.show()

exit()

fig, ax = plt.subplots(1)
plt.bar(np.arange(len(alphas)), meanAuc)
index = list(range(len(alphas)))
plt.xticks(index, alphas)
plt.show()

fig, ax = plt.subplots(1)
plt.bar(np.arange(len(alphas)), conceptPairs)
index = list(range(len(alphas)))
plt.xticks(index, alphas)
plt.show()
