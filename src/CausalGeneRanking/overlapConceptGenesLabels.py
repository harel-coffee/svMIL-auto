

import numpy as np
import sys
import matplotlib.pyplot as plt

conceptSG = np.loadtxt(sys.argv[1], dtype="object")
concept2 = np.loadtxt(sys.argv[2], dtype="object")
conceptPP = np.loadtxt(sys.argv[3], dtype="object")

allLabelsIntersect = list(set(conceptSG[:,0]) & set(concept2[:,0]) & set(conceptPP[:,0]))
SG2Intersect = list(set(conceptSG[:,0]) & set(concept2[:,0]))
SGPPIntersect = list(set(conceptSG[:,0]) & set(conceptPP[:,0]))
PP2Intersect = list(set(concept2[:,0]) & set(conceptPP[:,0]))

print "Number of genes in all labels: ", len(allLabelsIntersect)
print "Number of genes shared by SG and 2: ", len(SG2Intersect)
print "Number of genes shared by SG an PP: ", len(SGPPIntersect)
print "Number of genes shared by PP and 2: ", len(PP2Intersect)

np.savetxt('sharedConceptGenes.txt', allLabelsIntersect, delimiter='\t', fmt='%s')
exit()

from matplotlib_venn import venn3, venn3_circles

#3 categories, genes in COSMIC, genes with SNVs, genes that are DEG

v = venn3(subsets=(len(set(conceptSG[:,0])),len(set(concept2[:,0])),len(SG2Intersect), len(set(conceptPP[:,0])),len(SGPPIntersect), len(PP2Intersect),len(allLabelsIntersect)),
		  set_labels=('Somatic vs germline', '> 2 patients', 'Per patient'))
v.get_label_by_id('100').set_text(len(conceptSG[:,0]))
v.get_label_by_id('010').set_text(len(concept2[:,0]))
v.get_label_by_id('001').set_text(len(conceptPP[:,0]))
plt.title("")
plt.show()