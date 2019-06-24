"""
	Make a whoosh index with all Transcription factors per TAD
	Later, perhaps it is a good idea to also have other data types in here if that is faster than the current approach. 
	
	Documents are structured as:
	tad: chr_start_end
	tp: type of the element, e.g. tf for transcription factor
	start: start position of the element
	end: end position of the element
	
"""

# from whoosh.index import create_in
# from whoosh.fields import *
import os
import numpy as np
import sys
# from whoosh.query import *
# from whoosh.qparser import *
# 
# #0. Prepare the whoosh index to store data in
# 
# schema = Schema(tad=TEXT(stored=True), tp=TEXT(stored=True), start=NUMERIC(stored=True), end=NUMERIC(stored=True))
# if not os.path.exists("index"):
#     os.mkdir("index")
# ix = create_in("index", schema)
# writer = ix.writer()

#1. Read the TADs
tadFile = sys.argv[2]
tads = []
with open(tadFile, 'r') as inF:
	
	for line in inF:
		line = line.strip()
		splitLine = line.split("\t")
		
		tads.append([splitLine[0], int(splitLine[1]), int(splitLine[2])])

	tads = np.array(tads, dtype="object")

#2. Read the transcription factor file line by line and map it to the correct TAD
tfFile = sys.argv[1]
tadMap = dict()
lineCount = 0
with open(tfFile, 'r') as inF:
	for line in inF:
		
		if lineCount % 500000 == 0:
			print lineCount
		
		line = line.strip()
		
		splitLine = line.split("\t")
		
		tadChrSubset = tads[tads[:,0] == splitLine[0]]
		
		matchingStart = tadChrSubset[:,1] <= int(splitLine[1])
		matchingEnd = tadChrSubset[:,2] >= int(splitLine[2])
		
		matchingPos = matchingStart * matchingEnd
		if len(tadChrSubset[matchingPos]) < 1:
			continue
		matchingTad = tadChrSubset[matchingPos][0] #Take the first one for now if a tf is overlapping 2 TADs, celltypes are not the same between TADs and TFs
		
		tadName = matchingTad[0] + "_" + str(matchingTad[1]) + "_" + str(matchingTad[2])
		
		if tadName not in tadMap:
			tadMap[tadName] = []
		tadMap[tadName].append(splitLine[1])
		lineCount += 1
		
		#writer.add_document(tad=unicode(tadName), tp=u"tf", start=splitLine[1], end=splitLine[2])

outFile = sys.argv[3]
with open(outFile, 'w') as outF:
	for tad in tadMap:
		
		line = tad
		line += "\t"
		for tf in tadMap[tad]:
			line += str(tf)
			line += ","
		line += "\n"
		outF.write(line)

exit()		
# #Do a test search
# with ix.searcher() as searcher:
# 	parser=QueryParser("tad",schema=ix.schema)
# 	parser.add_plugin(GtLtPlugin())
# 	q=parser.parse(u"start:>40000 end:<500000")
# 	
# 	results = searcher.search(q)
# 	
# 	for result in results:
# 		print result
# 	
# exit()


# 
# writer.add_document(title=u"First document", path=u"/a", content=u"This is the first document we've added!")
# writer.add_document(title=u"Second document", path=u"/b", content=u"The second one is even more interesting!")

writer.add_document(tad=u"chr1_1_10000", tp=u"tf", chrom=u"chr1", start=50, end=55)
writer.add_document(tad=u"chr1_1_10000", tp=u"tf", chrom=u"chr1", start=70, end=75)

writer.commit()


with ix.searcher() as searcher:
	#query = And([Term("tad", u"chr1_1_10000"), Term("tf", u"chr1_50_55")])
	#query = QueryParser("tad", ix.schema).parse("chr1_1_10000")
	#query = And([Term("tad", u"chr1_1_10000"), Term("type", u"tf"), Term("chrom", u"1"), Term("start", ">30"), Term("end", "<60")])
	#query = And([Term("tad", u"chr1_1_10000"), Term("type", u"tf")])
	query = And([Term("tad", u"chr1_1_10000"), Term("start", 40), Term("end", 55), Term("chrom", u"chr1"), Term("tp", u"tf")])
	#query = And([Term("tad", u"chr1_1_10000"), Term("tf", u"chr1_50_55")])
	#results = searcher.search(query)
	
	
	parser=QueryParser("tad",schema=ix.schema)
	parser.add_plugin(GtLtPlugin())
	q=parser.parse(u"start:>40 end:<80 tad:chr1_1_10000")
	
	results = searcher.search(q)
	
	for result in results:
		print result
		
	
# from datetime import datetime, timedelta
# from whoosh import fields, index
# 
# schema = fields.Schema(title=fields.TEXT, content=fields.TEXT,
#                        date=fields.DATETIME)
# ix = index.create_in("index", schema)
# 
# w = ix.writer()
# w.add_document(title=u"Document 1", content=u"Rendering images from the command line",
#                date=datetime.utcnow())
# w.add_document(title=u"Document 2", content=u"Creating shaders using a node network",
#                date=datetime.utcnow() + timedelta(days=1))
# w.commit()
# 
# 
# from whoosh import index
# from whoosh.qparser import QueryParser
# from whoosh.qparser.dateparse import DateParserPlugin
# 
# ix = index.open_dir("index")
# 
# # Instatiate a query parser
# qp = QueryParser("content", ix.schema)
# 
# # Add the DateParserPlugin to the parser
# qp.add_plugin(DateParserPlugin())
# 
# qp = QueryParser("content", schema=ix.schema)
# 
# # Find all datetimes in 2005
# q = qp.parse(u"date:2005")
# 
# # Find all datetimes on June 24, 2005
# q = qp.parse(u"date:20050624")
# 
# # Find all datetimes from 1am-2am on June 24, 2005
# q = qp.parse(u"date:2005062401")
# 
# # Find all datetimes from Jan 1, 2005 to June 2, 2010
# q = qp.parse(u"date:[20050101 to 20100602]")
# 
# with ix.searcher() as searcher:
# 
# 	results = searcher.search(q)
# 	print results

