import matplotlib.pyplot as plt

outDir = sys.argv[1]

finalOutDir = outDir + '/methodComparison/results/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

#make a plot showing the true positive and false positive rates of each method.
#each sv type will get its own icon, and methods can be labeled by color

methods = ['chrCV MIL', 'chrCV simple RF', 'VEP', 'SVScore']
methodColors = ['#0055d4ff', '#c83737ff', '#9955ffff', '#808080ff']

tprsDEL = [0.64, 0.41, 0.02, 0.009]
fprsDEL = [0.19, 0.28, 0.2, 0.09]

tprsDUP = [0.52, 0.4, 0.08, 0.03]
fprsDUP = [0.33, 0.4, 0.46, 0.08]

tprsINV = [0.59, 0.27, 0, 0.007]
fprsINV = [0.24, 0.27, 0, 0.03]

tprsITX = [0.51, 0.44, 0, 0]
fprsITX = [0.24, 0.42, 0, 0.02]

#make the scatter plot
plt.scatter(fprsDEL, tprsDEL, marker='.', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsDUP, tprsDUP, marker='s', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsINV, tprsINV, marker='^', facecolor=methodColors, edgecolor=methodColors)
plt.scatter(fprsITX, tprsITX, marker='*', facecolor=methodColors, edgecolor=methodColors)

plt.savefig(finalOutDir + '/tpr_fpr.svg')



