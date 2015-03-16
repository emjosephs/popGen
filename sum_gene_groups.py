import sys

if len(sys.argv) < 3:
	print("python sum_gen_groups.py [data file] [gene list]")
	sys.exit()

# read in the gene list
geneList = []
genes = open(sys.argv[2],'r')
for gene in genes:
	geneList.append(gene.rstrip())

genes.close()

myData = open(sys.argv[1],'r')
header = myData.readline() #skip header
outData = [0]*(len(header.split()) - 1)

print(outData)

#for line in myData:
#	gene, data = line.split()[0], line.split()[1:]
#	if gene

