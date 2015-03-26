import sys
import summary
import random
import getopt
if len(sys.argv) < 4:
	print('python popgen_master.py [summary file] [output file] [pi/sfs/dfealpha]')
	sys.exit()
_d=320
def __main__():

	mySum = summary.Reader(open(sys.argv[1],'rb'))
	out = open(sys.argv[2],'w')
	#out.write("pac	piS	piN	Ssite#	Nsite#\n")
	myTest = sys.argv[3]
	processArgs(4)	
	count = _d
	siteDic = mySum.summary.typesAsInt()
	if myTest == "sfs":
		sfs(mySum, out, count, siteDic)

	elif myTest == "dfealpha":
		dfealpha(mySum, out, count, siteDic)
	elif myTest == "pi":
		dict = sfs(mySum, out, count, siteDic)
		pi_dict = pi(dict)
		print(pi_dict)
def dfealpha(mySum, out, count, siteDic):
	header = "gene	fold0.div	"+"	".join(["fold0."+str(x) for x in range(0,count)])+"	fold4.div	"+"	".join(["fold4."+str(x) for x in range(0,count)])
	sfsDict, divDict = {},{}
	for site in mySum:
		siteType = siteDic[int(site.Types[0])]
		#is there enough coverage?
		if site.TOTAL < _d:
			continue
		# is it the right type of site?
		if siteType not in ['0fold','4fold']:
			continue	
		# new gene?
		if site.GENE not in sfsDict.keys():
			sfsDict[site.GENE] = newSfsDict(count)
			divDict[site.GENE] = {'0fold':0, '4fold':0}

		siteMaf = maf(site.ALT_NUM, count)
		sfsDict[site.GENE][siteType][siteMaf] += 1 #add in sfs info
		#print(site.POS,sfsDict)
		#add in divergence info
		divDict[site.GENE][siteType] += site.DIVERGENCE

	#write out
	out.write(header)
	for gene in sfsDict.keys():
		out.write("\n"+gene)
		for thing in ['0fold','4fold']:
			out.write("	"+str(divDict[gene][thing]) + "	" + "	".join([str(x) for x in sfsDict[gene][thing]]))		

def sfs(mySum, out, count, siteDic):
	header = "gene	"+"	".join(["fold0."+str(x) for x in range(0,count)])+"	".join(["fold4."+str(x) for x in range(0,count)])
	geneDict = {}
	for site in mySum:
		siteType = siteDic[int(site.Types[0])]
		#is there enough coverage?
		if site.TOTAL < _d:
			continue
		# is it the right type of site?
		if siteType not in ['0fold','4fold']:
			continue	
		# new gene?
		if site.GENE not in geneDict.keys():
			#make empty dictionary 
			geneDict[site.GENE] = newSfsDict(count)
		
		#downsample and get alternate allele frequency
		#aaf = downsamp(site.REF_NUM, site.ALT_NUM,count)
		aaf = site.ALT_NUM
		#get minor allele count
		mafDown = maf(aaf, count)
		#add in this site info
		geneDict[site.GENE][siteType][mafDown] += 1
		#print(site.POS,site.GENE,siteType,mafDown,geneDict[site.GENE])	
	#write out
	out.write(header)
	for gene in geneDict.keys():
		out.write("\n"+gene)
		for thing in ['0fold','4fold']:
			out.write("	"+ "	".join([str(x) for x in geneDict[gene][thing]]))
	return(geneDict)
def pi (geneDict) :
	print(geneDict.keys())
	pi_gene = {}
	for gene in geneDict.keys():
		pi_gene[gene] =  {"0fold":0,"4fold":0}
		for frq in range(1,_d):
			pi_gene[gene]['0fold'] += (2*frq/(_d-1)*(1-frq/(_d-1)))*geneDict[gene]['0fold'][frq]
			pi_gene[gene]['4fold'] += (2*frq/(_d-1)*(1-frq/(_d-1)))*geneDict[gene]['4fold'][frq]
	return(pi_gene)
def downsamp(ref, alt, count): #takes in counts of ref and alt alleles, and count=N of downsampling, returns the # of alt alleles in the downsampled sample
	totList = [0]*ref + [1]*alt
	ds = random.sample(totList, count)
	return(sum(ds))


def newSfsDict(depth):
	mySfs = {"0fold":[],"4fold":[]}
	for i in range(0,depth):
		mySfs["0fold"].append(0)
		mySfs["4fold"].append(0)
	return(mySfs)	


def maf(aaf, count): # calculates the # of minor alleles when given the # of alternate alleles and the total sample count
	if aaf > float(count)/2:
		maf = count - aaf
	else:
		maf = aaf
	return(maf)


def processArgs(num):
	try:
		opts, args = getopt.getopt(sys.argv[num:],"d:")
	except getopt.GetoptError:
		"wrong usage"
	for opt,arg in opts:
		# depth cutoff
		if opt=="-d":   
			global _d
			_d=int(arg)
#def theta_w (
if __name__ == "__main__":
        __main__()





