import pandas as pd
import patsy
import sys
import numpy as np
from combat import combat

#Read the normalized rna-seq data

normData = pd.read_csv(sys.argv[1], index_col=0, sep='\t')

print(normData)

#make a table for phenotype data

#batch labels
batches = ['GTEx']*479 + ['BRCA']*162
tissue = ['Normal']*479 + ['Tumor']*162

pheno = [batches, tissue]
pheno = np.array(pheno, dtype='object')
phenoData = pd.DataFrame(pheno).T
phenoData.rows = list(normData.columns.values)
phenoData.columns = ['Batch', 'Tissue']

print(phenoData)

#make the 'mod' thing, whatever it may be
mod = patsy.dmatrix("~ Tissue", phenoData, return_type="dataframe")


ebat = combat(normData, batches, mod)

print(ebat)


