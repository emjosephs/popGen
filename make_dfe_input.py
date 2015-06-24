import sys
import random

if len(sys.argv) < 4:
	print('python make_dfe_input.py [input prefix] [gene list] [label] [boot?]') 
	sys.exit()
if len(sys.argv) < 5:
	dowebootstrap = False
elif sys.argv[4] == "boot":
	dowebootstrap = True

def __main__():

	#make list of gene names
	geneList = [] 
	genes = open(sys.argv[2],'r')
	for gene in genes:
        	geneList.append(gene.rstrip())
	genes.close()

	if dowebootstrap: #we bootstrap
		bootList = geneList
		geneList = bootstrap(bootList)

	modName = sys.argv[3]
	numFields = None

	for i in range(1,9):
        	myData = open(sys.argv[1]+str(i),'r')
        	header = myData.readline() #skip header
        	if not numFields: #if first file, figure out the number of fields 
			outData, dataNames = makeOutData(header)
			numFields = len(dataNames)		
	
        	for line in myData:
                	if line[0:4] == "gene": #skip errant headers
                        	continue 
                	gene, data = line.split()[0], [int(x) for x in line.split()[1:]]
                	if gene in geneList: #is it in our module
                        	for i in range(0, numFields):
                                	outData[i].append(data[i])

	outDic = dict( zip( dataNames,outData) )	
	selSFSNames = ["fold0."+str(x) for x in range(0, 319)]
	neuSFSNames = ["fold4."+str(x) for x in range(0, 319)]
	selSFS = [sum(outDic[x]) for x in selSFSNames]
	neuSFS = [sum(outDic[x]) for x in neuSFSNames]
	selDiv, neuDiv = sum(outDic["fold0.div"]), sum(outDic["fold4.div"])

	#print out in dfe alpha input format
	print(modName) #[NAME1]
	print(str(sum(selSFS)) +"	"+ str(selDiv)) #[SEL SITES] [SEL DIFFS]
	print( str(sum(neuSFS)) +"	" + str(neuDiv))  #[NEU SITES] [NEU DIFFS]
	print("320") #[d1]
	print( "	".join([str(x) for x in selSFS]) + "	0")#[Selected SFS 1]
	print( "	".join([str(x) for x in neuSFS]) +"	0")#[Neutral SFS 1]


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

def bootstrap(geneList): #does a bootstrap by sampling with replacement
	outList = []
	for i in range(0,len(geneList)):
		outList.append(random.choice(geneList))
	return(outList)

if __name__ == "__main__":
        __main__()



