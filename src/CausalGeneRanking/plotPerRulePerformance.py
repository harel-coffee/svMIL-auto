"""
	Get the scores for every feature, and plot these with vertical lines in a plot

"""

import matplotlib.pyplot as plt
import numpy as np

#Gains

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [5.96768130283696e-75, 9.025675391175049e-72, 0.0698324677480698, 0.07607505848668761, 8.073202703436674e-43, 2.1618844504538086e-33, 0.015465298855195115, 5.264953341705389e-44]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_cosmic.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [1.7326360323860752e-92, 0.0, 0.9689763166970913, 0.9890939248758223, 0.009392785206533907, 0.0022061247071801537, 0.8513044985435061, 0.011692728479067016]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_bc.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0.0, 0.0, 8.616640136721214e-06, 8.183109137078339e-05, 7.293214939150255e-157, 4.682793945369012e-117, 1.3973022059287285e-07, 4.5088936248048134e-144]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg.svg')

#Repeat for lost features

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [9.528605281531696e-13, 0.05978811082646416, 0.02570014473785348, 0.08221628121775268, 2.1602575486244998e-54, 0.03954706422738248, 0.003802621429074631, 1.7513477811461663e-67]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_cosmic.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0.0009987757824814457, 1, 0.3029360733233496,  0.293353510183918, 0.0020952403217837415, 0.3846460321577463, 0.3316850680568625, 1.2088355114059256e-06]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_bc.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [1.2596309965999392e-41, 3.238976906145921e-07, 4.209544493298611e-05, 6.237975326899646e-08, 2.486408621376771e-171, 5.064175661770919e-08, 7.571272442603699e-07, 3.3154844222239676e-177]
y = [i + 1e-250 for i in y]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, -np.log(y))
plt.axhline(y=-np.log(0.05), color='k')
plt.title('P-values for genes that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('-log(p-value)')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg.svg')

#repeat but then make the plot for the total number of genes (counts) and add the average counts in shuffled case


#Gains

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [6, 5, 7, 8, 11, 11, 11, 11]
x = range(0, len(y))

fig, ax = plt.subplots()
ax.stem(x, y)
averages = [28, 29, 4, 5, 53, 53, 7, 52]
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+1)
plt.title('Genes that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_cosmic_counts.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [1, 1, 1, 1, 1, 1, 1, 1]
x = range(0, len(y))
averages = [3, 3, 1, 1, 3, 3, 1, 3]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+1)
plt.title('Genes that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_bc_counts.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [37, 32, 64, 65, 87, 87, 86, 87]
x = range(0, len(y))
averages = [562, 561, 41, 45, 605, 602, 57, 598]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+10)
plt.title('Genes that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Gained features')
plt.tight_layout()
plt.savefig('gains_deg_counts.svg')

#losses

#First for COSMIC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0, 5, 3, 10, 15, 15, 15, 0]
x = range(0, len(y))

averages = [12, 2, 6, 7, 73, 10, 10, 75]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+2)
plt.title('Genes that are DEG and in COSMIC')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_cosmic_counts.svg')
#plt.show()

#Repeat for BC + DEGs

#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [0, 0, 0, 0, 2, 2, 2, 0]
x = range(0, len(y))
averages = [2, 0, 1, 1, 5, 1, 1, 5]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+1)
plt.title('Genes that are DEG and breast cancer genes')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_bc_counts.svg')

#For DEGs only
#eqtls, enhancers, promoters, CpG, TFs, HiC, histones, dnaseI
y = [47, 60, 83, 91, 111, 110, 108, 111]
x = range(0, len(y))
averages = [187, 38, 58, 62, 780, 77, 76, 781]

fig, ax = plt.subplots()
ax.stem(x, y)
count = 0
for i,j in zip(x,y):
	ax.annotate(str(averages[count]), xy=(i,j), xytext=(-2,10), textcoords='offset points') #add the averages
	count += 1
plt.ylim(np.min(y), np.max(y)+10)
plt.title('Genes that are DEG')
ax.set_xticklabels(['', 'eQTLs', 'enhancers', 'promoters', 'CpG', 'TFs', 'HiC', 'Histones', 'DNAseI'])
plt.xticks(rotation=90)
plt.ylabel('Number of genes')
plt.xlabel('Lost features')
plt.tight_layout()
plt.savefig('losses_deg_counts.svg')
