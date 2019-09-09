"""
	Get the scores for every feature, and plot these with vertical lines in a plot

"""

import matplotlib.pyplot as plt
import numpy as np

#Gains

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0.0, 0.0, 1.74e-20, 6.02e-20, 0.0, 0.0, 3.15e-37, 0.0]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_cosmic.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [2.37e-200, 2.36e-200, 1.18e-05, 6.92e-05, 4.66e-55, 4.46e-54, 1.09e-08, 1.87e-42]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_bc.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0.0, 0.0, 3.78e-74, 6.97e-65, 0.0, 0.0, 1.79e-83, 0.0]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg.svg')

#Gains

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [29, 27, 8, 9, 68, 68, 13, 67]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, y)
averages = [0, 0, 1, 1, 2, 2, 1, 1]
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+10)
plt.title('SV-gene pairs that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_cosmic_counts.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [3, 3, 1, 1, 5, 5, 1, 5]
x = range(0, len(y))
averages = [0, 0, 0, 0, 0, 0, 0, 0]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+1)
plt.title('SV-gene pairs that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_bc_counts.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [578, 579, 74, 76, 709, 705, 96, 705]
x = range(0, len(y))
averages = [5, 5, 12, 12, 32, 29, 17, 29]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+50)
plt.title('SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_counts.svg')


#Repeat for lost features

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [5.07e-223, 7.16e-34, 1.04e-40, 1.78e-25, 0.0, 3.74e-45, 2.77e-69, 0.0]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_cosmic.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [2.37e-200, 1, 0.82, 0.89, 3.11e-77, 3.85e-31, 2.92e-14, 1.31e-57]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_bc.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0.0, 4.00e-82, 3.63e-83, 2.04e-87, 0.0, 9.7e-116, 4.83e-87, 0.0]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg.svg')

#repeat but then make the plot for the total number of genes (counts) and add the average counts in shuffled case




#losses

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [18, 9, 12, 12, 97, 18, 18, 98]
x = range(0, len(y))

averages = [0, 1, 1, 1, 2, 1, 1, 2]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+10)
plt.title('SV-gene pairs that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_cosmic_counts.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [3, 0, 0, 0, 6, 2, 2, 6]
x = range(0, len(y))
averages = [0, 0, 0, 0, 0, 0, 0, 0]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+1)
plt.title('SV-gene pairs that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_bc_counts.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [229, 71, 88, 97, 938, 126, 124, 942]
x = range(0, len(y))
averages = [11, 10, 16, 17, 41, 22, 22, 41]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+100)
plt.title('SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_counts.svg')

#plot combination of features

#losses

#DEGs
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [229, 93, 105, 115, 127, 127, 127, 127]
x = range(0, len(y))
averages = [X, 615, 729, 3200, 4105, 4118, 4118, 4118] #total number of pairs

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+100)
plt.ylim(np.min(y), np.max(y))
plt.title('SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of SV-gene pairs')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_counts_addition.svg')

#gains

#DEGs
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [578, 54, 79, 83, 97, 97, 97, 97]
x = range(0, len(y))
averages = [X, 1525, 2377, 2592, 3174, 3174, 3174, 3175] #total number of pairs

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+100)
plt.ylim(np.min(y), np.max(y))
plt.title('SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of SV-gene pairs')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_counts_addition.svg')

#gains + losses

#DEGs
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [70, 97, 110, 118, 127, 127, 127, 127]
x = range(0, len(y))
averages = [2433, 3280, 3949, 4251, 4692, 4693, 4693, 4708] #total number of pairs

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+100)
plt.ylim(np.min(y), np.max(y))
plt.title('SV-gene pairs that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of SV-gene pairs')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_counts_addition.svg')
