import sys
import summary
import random
import getopt
import math
if len(sys.argv) < 4:
	print('python popgen_master.py [summary file] [output file] [diversity/sfs/dfealpha]')
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
	elif myTest == "diversity":
		diversity(mySum, out, count, siteDic)
	


def dfealpha(mySum, out, count, siteDic):
	header = "gene	fold0.div	"+"	".join(["fold0."+str(x) for x in range(0,count)])+"	fold4.div	"+"	".join(["fold4."+str(x) for x in range(0,count)])
	sfsDict, divDict = {},{}
	out.write(header)
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
def diversity (mySum, out, count, siteDic):
	divDict = {}
	out.write("gene\tn_rep\tpoly_rep\ttotal_rep\ttheta_rep\tpi_rep\tTajimas_D_rep\tn_syn\tpoly_syn\ttotal_syn\ttheta_syn\tpi_syn\tTajimas_D_syn\n")
	for site in mySum:
	        siteType = siteDic[int(site.Types[0])]
                #is there enough coverage?
	        if site.TOTAL < _d:
        	 	  continue
        	        # is it the right type of site?
        	if siteType not in ['0fold','4fold']:
        	          continue
                # new gene?
        	if site.GENE not in divDict.keys():
                        #make empty dictionary
        	          divDict[site.GENE] = {'0fold':{'poly':0,'total':0,'pi':0,'n':0},'4fold':{'poly':0,'total':0,'pi':0,'n':0}}
		if site.REF_NUM < site.TOTAL and site.ALT_NUM < site.TOTAL:
			divDict[site.GENE][siteType]['poly']+=1
		divDict[site.GENE][siteType]['total']+=1
		aaf = float(site.ALT_NUM)
		n = float(site.TOTAL)
		divDict[site.GENE][siteType]['n']+= site.TOTAL
		divDict[site.GENE][siteType]['pi']+=float(2*aaf/n)*(1-aaf/n)*(n/(n-1)) # see Gillespie 2nd ed p. 45
#	print divDict
	genelist = divDict.keys()
	genelist.sort()
	for gene in  genelist:
		out.write(gene+"\t")
		for type in ['0fold','4fold']:
			theta = ''
			a = ''
			pi = ''
			tD=''
			if divDict[gene][type]['total']>0:
				mean_n=int(divDict[gene][type]['n']/divDict[gene][type]['total'])
				a = calc_a(mean_n)  #Watterson's thera : S/a
		#		print(gene,type,a,mean_n)
				theta=(divDict[gene][type]['poly']/a)/divDict[gene][type]['total']
				pi=divDict[gene][type]['pi']/divDict[gene][type]['total']
		#		print (a,theta,pi)
				den = D_denominator(mean_n,divDict[gene][type]['poly'])
		#		print("den_syn",den)
				if den > 0:
					tD = (pi-theta)/den
#					print (gene, den, tD)
				sep="\t"
			out.write(sep.join(str(x) for x in [mean_n,divDict[gene][type]['poly'],divDict[gene][type]['total'],theta,pi,tD]))
			out.write("\t")
		out.write("\n")
def calc_a(n):
	a = float(0)
	for i in range(1,n):
		a+= float(1.0/i)
#		print (i,n)
	return(a)
def D_denominator(N,S): # denominator for Tajima's D
   a1 =0.0 #denom of theta
   for i in range (1,N):
      a1=a1 + 1.0/i
   

   a2 = 0.0 #denom of theta
   for i in range (1,N):
       a2=a2 + 1.0/(i*i)
    

   b1=(N+1.0)/(3.0*(N-1.0))
   b2=(2.0*(N*N+N+3.0))/(9.0*N*(N-1.0))
   c1=b1-(1.0/a1)
   e1=c1/a1
   c2 = b2-(N+2.0)/(a1*N) + a2/(a1*a1)
   e2 = c2/((a1*a1)+a2)
   S= float(S)
   denominator = math.sqrt(e1*S+e2*S*(S-1.0))
   return(denominator)
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





