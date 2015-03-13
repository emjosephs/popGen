# adds gene names to summary file
import summary
import sys

gffDic = {} #keys are scaffolds, values are a list of positions, 0 for 0s, or else gene name
for i in range(1,9):
	gffDic[i] = [0]*19624520
#	gffDic[i] = [0]*2000

#make a dic of all the gene locations
geneFile = open('../gff_gene_flatfile','r')
#geneFile = open('../test_gff_gene_flatfile','r')
for line in geneFile:
	scaf, gene, start, end = line.split()
	scafNum = int(scaf.split('_')[1])
	for i in range(int(start), int(end)+1): #not that worried about 0 or 1 based coordinates because we're only dealing with 0 and 4 fold sites.
		gffDic[scafNum][i] = gene

geneFile.close()

#myOut = open('../../0fold4fold.withgenes.summary','w'))
sumRead = summary.Reader(open('../../0fold4fold.summary','rb'))
#sumRead = summary.Reader(open('../test.summary','rb'))

sumRead.addGenes()
print(sumRead.Header)
for site in sumRead:
	if gffDic[int(site.CHROM.split('_')[1])][site.POS] == 0: #not in a gene?
		sys.stderr.write("err not in gene "+site.prettyStr()+"\n")
		continue
	else:
		geneName = gffDic[int(site.CHROM.split('_')[1])][site.POS]
		print("	".join([site.__str__(),geneName,"NA"]))
	




