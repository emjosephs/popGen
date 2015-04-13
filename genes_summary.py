# adds gene names to summary file
import summary
import sys

if len(sys.argv) < 3:
	print('python gene_summary.py [input summary file]')
	sys.exit()

gffDic = {} #keys are scaffolds, values are a list of positions, 0 for 0s, or else gene name
for i in range(1,9):
	gffDic[i] = [0]*19624520
#	gffDic[i] = [0]*2000

#make a dic of all the gene locations
geneFile = open(sys.argv[2],'r')
#geneFile = open('../test_gff_gene_flatfile','r')
for line in geneFile:
	scaf, gene, start, end = line.split()
	scafNum = int(scaf.split('_')[1])
	for i in range(int(start), int(end)+1): #not that worried about 0 or 1 based coordinates because we're only dealing with 0 and 4 fold sites.
		gffDic[scafNum][i] = gene

geneFile.close()

#myOut = open('../../0fold4fold.withgenes.summary','w'))
#sumRead = summary.Reader(open('/data/youngwha.lee/189_genomes/UG_all_vars/recal_vcfs/vcfsummaries/downsampled320/sc8_down_320_4','rb'))
sumRead = summary.Reader(open(sys.argv[1],'rb'))

sumRead.addGenes()
print(sumRead.Header)
for site in sumRead:
	if gffDic[int(site.CHROM.split('_')[1])][site.POS] == 0: #not in a gene?
		sys.stderr.write("err not in gene "+site.prettyStr()+"\n")
		continue
	else:
		geneName = gffDic[int(site.CHROM.split('_')[1])][site.POS]
		site.GENE = geneName
		print(site)
		#print("	".join([site.__str__(),geneName,"NA"]))
	




