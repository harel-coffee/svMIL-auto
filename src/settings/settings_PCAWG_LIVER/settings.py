files = dict(

	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty.
	svFile = '../../data/svs/liver_pcawg_parsed.txt',
	causalGenesFile = '../../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	promoterFile = '../../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../../data/tf/tf_experimentallyValidated.bed', #All validated human TFs
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../../data/chromosomes/hg19Coordinates.txt",
	snvDir = '../../data/snvs/snv_mnv/',
	cnvFile = '../../data/cnvs/all_samples.consensus_CN.by_gene.170214.txt',
	geneNameConversionFile = '../../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs, converting ENSG identifiers to gene names.
	metaDataFile = '../../data/svs/icgc_metadata.tsv', #metadata to map identifiers
	expressionFile = '../../data/expression/tophat_star_fpkm_uq.v2_aliquot_gl.tsv',

	#specific for liver, mixed with HMEC where data unavailable
	eQTLFile = '../../data/eQTLs/liver/liver_eQTLs.bed', #File with eQTLs specific for liver tissue
	tadFile = "../../data/tads/liver/Liver_STL011_Leung2015-raw_TADs.bed", #File with TADs specific for liver tissue
	enhancerFile = '../../data/enhancers/liver/encoderoadmap_elasticnet.64.csv',
	h3k9me3 = '../../data/histones/liver/ENCFF845SNW_H3K9me3.bed',
	h3k4me3 = '../../data/histones/liver/ENCFF065NRN_H3K4me3.bed',
	h3k27ac = '../../data/histones/liver/ENCFF752QSK_H3K27ac.bed',
	h3k27me3 = '../../data/histones/liver/ENCFF933NEK_H3K27me3.bed',
	h3k4me1 = '../../data/histones/liver/ENCFF327BBS_H3K4me1.bed',
	h3k36me3 = '../../data/histones/liver/ENCFF744HEP_H3K36me3.bed',
	dnaseIFile = '../../data/dnase/liver/ENCFF286LYP.bed',
	chromHmmFile = '../../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed', #HMEC
	rnaPolFile = '../../data/rnapol/hmec/ENCFF433ZKP.bed', #HMEC
	superEnhancerFile = '../../data/superEnhancers/liver/se_20200411_hepg2.bed',
	ctcfFile = '../../data/ctcf/liver/ENCFF690BYG.bed',
	hicFile = '../../data/hic/HMEC_groupedTADInteractions.txt', #HMEC

)

general = dict(

	source = 'PCAWG',
	cancerType = 'LIVER',
	shuffleTads = False #Should TAD positions be shuffled
)
