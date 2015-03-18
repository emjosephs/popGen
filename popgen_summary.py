import sys
import summary
import random

if len(sys.argv) < 4:
	print('python popgen_master.py [summary file] [output file] [pi/sfs/dfealpha]')
	sys.exit()

def __main__():

	mySum = summary.Reader(open(sys.argv[1],'rb'))
	out = open(sys.argv[2],'w')
	#out.write("pac	piS	piN	Ssite#	Nsite#\n")
	myTest = sys.argv[3]
	count = 320
		
	siteDic = mySum.summary.typesAsInt()
	if myTest == "sfs":
		sfs(mySum, out, count, siteDic)

	elif myTest == "dfealpha":
		dfealpha(mySum, out, count, siteDic)

def dfealpha(mySum, out, count, siteDic):
	header = "gene	fold0.div	"+"	".join(["fold0."+str(x) for x in range(0,count)])+"	fold4.div	"+"	".join(["fold4."+str(x) for x in range(0,count)])
	sfsDict, divDict = {},{}
	for site in mySum:
		siteType = siteDic[int(site.Types[0])]
		#is there enough coverage?
		if site.TOTAL < 320:
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

		#add in divergence info
		divDict[site.GENE][siteType] += site.DIVERGENCE

	#write out
	out.write(header)
	for gene in sfsDict.keys():
		out.write("\n"+gene)
		for thing in ['0fold','4fold']:
			out.write("	"+str(divDict[thing]) + "	" + "	".join([str(x) for x in sfsDict[gene][thing]]))		

def sfs(mySum, out, count, siteDic):
	header = "gene	"+"	".join(["fold0."+str(x) for x in range(0,count)])+"	".join(["fold4."+str(x) for x in range(0,count)])
	for site in mySum:
		geneDict = {}
		siteType = siteDic[int(site.Types[0])]
		#is there enough coverage?
		if site.TOTAL < 320:
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
			
	#write out
	out.write(header)
	for gene in geneDict.keys():
		out.write("\n"+gene)
		for thing in ['0fold','4fold']:
			out.write("	"+ "	".join([str(x) for x in geneDict[gene][thing]]))


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


if __name__ == "__main__":
        __main__()





