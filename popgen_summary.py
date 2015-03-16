import sys
import summary
import random

if len(sys.argv) < 4:
	print('python popgen_master.py [summary file] [output file] [pi/sfs]')
	sys.exit()

def __main__():

	mySum = summary.Reader(open(sys.argv[1],'rb'))
	out = open(sys.argv[2],'w')
	#out.write("pac	piS	piN	Ssite#	Nsite#\n")
	myTest = sys.argv[3]
	geneDict = {}
	count = 160

		
	if myTest == "sfs":
		siteDic = mySum.summary.typesAsInt()
		header = "gene	"+"	".join(["0fold."+str(x) for x in range(0,count)])+"	".join(["4fold."+str(x) for x in range(0,count)])
		for site in mySum:
			siteType = siteDic[int(site.Types[0])]
			#is there enough coverage?
			if site.TOTAL < 161:
				continue
			# new gene?
			if site.GENE not in geneDict.keys():
				#make empty dictionary 
				geneDict[site.GENE] = newSfsDict(count)
		
			#downsample and get alternate allele frequency
			aaf = downsamp(site.REF_NUM, site.ALT_NUM,count)
			#get minor allele count
			mafDown = maf(aaf, count)
			#add in this site info
			geneDict[site.GENE][siteType][mafDown] += 1
			
		#write out
		out.write(header)
		for gene in geneDict.keys():
			out.write("\n"+gene+"	")
			for thing in ['0fold','4fold']:
				out.write( "	".join([str(x) for x in geneDict[gene][thing]]))


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





