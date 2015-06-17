#returns allele frequency at all sites as a function of distance from gene for each gene

import sys
import summary_distance
from scipy.stats import sem

if len(sys.argv) < 2:
	print('python gene_distance_maf.py [summary file]')
	sys.exit() 

def __main__():

	
	upList = [] 
	dwList = []
	for i in range(0,50):
		upList.append([]) #bins of 100 bp
		dwList.append([])

	# read through summary file
	mysum = summary_distance.Reader(open(sys.argv[1],'rb'))
	for record in mysum:
		#is it polymorphic?
		if record.ALT_NUM == 0:
			continue

		#calculate MAF
		myMaf = maf320(record.REF_NUM, record.ALT_NUM)
		
		if record.UPSTR[0] != "NA":#add to upstr list
			for gene in record.UPSTR:
				geneDist = int(gene.split(':')[1])
				binNo = binIt(geneDist)
				upList[binNo].append(myMaf)
		
		if record.DWSTR[0] != "NA": #add to dwstr list
			for gene in record.DWSTR:
				geneDist = int(gene.split(':')[1])
				binNo = binIt(geneDist)
				dwList[binNo].append(myMaf)

	#calc means
	upMeans = [float(sum(x))/len(x) for x in upList]
	dwMeans = [float(sum(x))/len(x) for x in dwList]
	upSte = [sem(x) for x in upList]	
	dwSte = [sem(x) for x in dwList]	


	#print out
	print("dir	"+"	".join([str(x) for x in range(0,50)]))
	print("upMeans	"+"	".join([str(x) for x in upMeans]))
	print("dwMeans	"+"	".join([str(x) for x in dwMeans]))
	print("upSte	"+"	".join([str(x) for x in upSte]))
	print("dwSte	"+"	".join([str(x) for x in dwSte]))


def maf(ref,alt):
	propRef = float(ref)/(ref+alt)
	if propRef <= 0.5:
		return(propRef)
	else:
		return(1-propRef)

def maf320(ref,alt):
	if ref >= 160:
		return(alt)
	else:
		return(ref)

def binIt(nu):
	if nu == 0:
		return(0)
	elif nu == 5000:
		return(49)
	else:
		return(int(abs(round((float(nu)-50)/100))))



if __name__ == "__main__":   
    __main__()
	
