Folder with genomic positions of Hi-C interactions

Preprocessing steps are also in preprocess.sh.

Sources:
HMEC_interactions.txt: file with all Hi-C interactions in the HMEC tissue. Processed from the 5kb bins from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63525&format=file&file=GSE63525%5FHMEC%5Fintrachromosomal%5Fcontact%5Fmatrices%2Etar%2Egz
using DataProcessing/Parsers/HiC/thresholdHiCData.py. The output of this script was copied to data/thesholdedHiC. Then ran DataProcessing/parseHicInteractions.py and ran DataProcessing/makeTADHiCFile.py on the output of this file.
HMEC_groupedTADInteractions.txt: processed from HMEC_interactions.txt and tads/HMEC_Lieberman-raw_TADs.bed using DataProcessing/makeTADHiCFile.py
