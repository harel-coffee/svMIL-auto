import matplotlib.pyplot as plt
import sys
import os

outDir = sys.argv[1]

finalOutDir = outDir + '/figure3d/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

#make a plot showing the true positive and false positive rates of each method.
#each sv type will get its own icon, and methods can be labeled by color

methods = ['chrCV MIL', 'chrCV simple RF', 'VEP', 'SVScore']
methodColors = ['#0055d4ff', '#c83737ff', 'orange', '#808080ff']

#These are obtained from running the individual scripts for each method (see workflow.sh)
tprsDEL = [0.53, 0.56, 0.02, 0.009]
fprsDEL = [0.20, 0.56, 0.2, 0.09]

tprsDUP = [0.58, 0.45, 0.08, 0.03]
fprsDUP = [0.30, 0.46, 0.46, 0.08]

tprsINV = [0.60, 0.38, 0, 0.007]
fprsINV = [0.25, 0.37, 0, 0.03]

tprsITX = [0.62, 0.47, 0, 0]
fprsITX = [0.30, 0.43, 0, 0.02]

#make the scatter plot
plt.scatter(fprsDEL, tprsDEL, marker='.', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsDUP, tprsDUP, marker='s', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsINV, tprsINV, marker='^', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsITX, tprsITX, marker='*', facecolor=methodColors, edgecolor=methodColors)

plt.savefig(finalOutDir + '/tpr_fpr.svg')



