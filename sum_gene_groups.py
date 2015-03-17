import sys
import numpy

if len(sys.argv) < 3:
	print("python sum_gen_groups.py [data file] [gene list]")
	sys.exit()

# read in the gene list
geneList = []
genes = open(sys.argv[2],'r')
for gene in genes:
	geneList.append(gene.rstrip())

genes.close()
modName = sys.argv[2].split('/')[-1].split('.')[0]

myData = open(sys.argv[1],'r')
header = myData.readline() #skip header
outData = []
numFields = len(header.split())-1
for i in range(0,numFields):
	outData.append([])
print(numFields)

for line in myData:
	gene, data = line.split()[0], [int(x) for x in line.split()[1:]]
	if gene in geneList:
		#outData = [x+y for x,y in zip(outData,data)]
		for i in range(0, numFields):
			outData[i].append(data[i])
means = [float(sum(x))/len(x) for x in outData]
print(modName+"	means	"+means)
errs =  [numpy.std(numpy.array(x)) for x in outData]
print(modName+"	stdevs	"+errs)



