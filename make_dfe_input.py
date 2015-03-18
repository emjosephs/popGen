import sys

if len(sys.argv) < 4:
	print('python make_dfe_input.py [input prefix] [gene list] [label]') 
	sys.exit()

def __main__():

	#make list of gene names
	geneList = [] 
	genes = open(sys.argv[2],'r')
	for gene in genes:
        	geneList.append(gene.rstrip())
	genes.close()

	modName = sys.argv[3]
	numFields = None

	for i in range(1,9):
        	myData = open(sys.argv[1]+str(i),'r')
        	header = myData.readline() #skip header
        	if not numFields: #if first file, figure out the number of fields 
			outData, dataNames = makeOutData(header)
			
        	for line in myData:
                	if line[0:4] == "gene": #skip errant headers
                        	continue 
                	gene, data = line.split()[0], [int(x) for x in line.split()[1:]]
                	if gene in geneList: #is it in our module
                        	for i in range(0, numFields):
                                	outData[i].append(data[i])

	#print out in dfe alpha input format
	


#means = [str(float(sum(x))/len(x)) for x in outData]
#print(modName+" means   "+"     ".join(means))
#errs =  [str(numpy.std(numpy.array(x))) for x in outData]
#print(modName+" stdevs  "+"     ".join(errs))

def makeOutData(header):
	numFields = len(header.split())-1
        outData = []
        for i in range(0,numFields):
		outData.append([]) #make empty outData
	dataNames = header.split()[1:]  #make a list of names for parsing the outdata later
	return(outData, dataNames)

if __name__ == "__main__":
        __main__()



