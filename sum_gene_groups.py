import sys
import numpy

if len(sys.argv) < 3:
	print("python sum_gen_groups.py [data file prefix] [gene list]")
	sys.exit()

# read in the gene list
geneList = []
genes = open(sys.argv[2],'r')
for gene in genes:
	geneList.append(gene.rstrip())

genes.close()
modName = sys.argv[2].split('/')[-1].split('.')[0]
numFields = None

for i in range(1,9):
	myData = open(sys.argv[1]+str(i),'r')
	header = myData.readline() #skip header
	if not numFields: #if first file, figure out the number of fields 
		numFields = len(header.split())-1
		outData = []
		for i in range(0,numFields):
			outData.append([]) #,ale pitData
	for line in myData:
		if line[0:4] == "gene":
			continue #header gile
		gene, data = line.split()[0], [int(x) for x in line.split()[1:]]
		if gene in geneList:
			#outData = [x+y for x,y in zip(outData,data)]
			for i in range(0, numFields):
				outData[i].append(data[i])

print(header)
means = [str(float(sum(x))/len(x)) for x in outData]
print(modName+"	means	"+"	".join(means))
errs =  [str(numpy.std(numpy.array(x))) for x in outData]
print(modName+"	stdevs	"+"	".join(errs))



