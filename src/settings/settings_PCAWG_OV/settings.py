files = dict(

	svFile = '../data/svs/ov_pcawg_parsed.txt',
	causalGenesFile = '../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	promoterFile = '../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../data/tf/tf_experimentallyValidated.bed', #All validated human TFs
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../data/chromosomes/hg19Coordinates.txt",
	snvDir = '../data/snvs/snv_mnv/',
	cnvFile = '../data/cnvs/all_samples.consensus_CN.by_gene.170214.txt',
	geneNameConversionFile = '../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs, converting ENSG identifiers to gene names.
	metaDataFile = '../data/svs/icgc_metadata.tsv', #metadata to map identifiers
	expressionFile = '../data/expression/tophat_star_fpkm_uq.v2_aliquot_gl.tsv',

	#specific for ov, mixed with HMEC where data unavailable
	eQTLFile = '../data/eQTLs/ovarian_eQTLs.bed', #File with eQTLs specific for ovarian tissue
	tadFile = "../data/tads/ov/ovary_hg38_liftover.bed", #File with TADs specific for ovarian tissue
	enhancerFile = '../data/enhancers/ov/encoderoadmap_elasticnet.95.csv',
	h3k9me3 = '../data/histones/ov/ENCFF717WXC_H3K9me3.bed',
	h3k4me3 = '../data/histones/ov/ENCFF320JHG_H3K4me3.bed',
	h3k27ac = '../data/histones/ov/ENCFF657AUA_H3K27ac.bed',
	h3k27me3 = '../data/histones/ov/ENCFF712UCB_H3K27me3.bed',
	h3k4me1 = '../data/histones/ov/ENCFF917PWI_H3K4me1.bed',
	h3k36me3 = '../data/histones/ov/ENCFF302DXB_H3K36me3.bed',
	dnaseIFile = '../data/dnase/ov/ENCFF883WWT.bed',
	chromHmmFile = '../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed', #HMEC
	rnaPolFile = '../data/rnapol/ov/ENCFF570SMG.bed',
	superEnhancerFile = '../data/superEnhancers/ov/se_20200408_ovary.bed',
	ctcfFile = '../data/ctcf/ov/ENCFF522DLJ.bed',
	hicFile = '../data/hic/HMEC_groupedTADInteractions.txt',

)

general = dict(

	source = 'PCAWG',
	cancerType = 'OV',
	shuffleTads = False #Should TAD positions be shuffled
)

