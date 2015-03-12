import sys
import vcf
import random

if len(sys.argv) < 4:
	print('python popgen_master.py [list of vcfs] [output file] [pi/sfs]')
	sys.exit()

def __main__():

	vcfList = open(sys.argv[1],'r')
	out = open(sys.argv[2],'w')
	#out.write("pac	piS	piN	Ssite#	Nsite#\n")
	myTest = sys.argv[3]

	for line in vcfList:
		myVcf = open(line.rstrip(),'r')
		pac = line.split('/')[3].split('.')[0]
		
		if myTest == "pi":
			piDict = {"syn":[], "non":[]}
			for record in vcf.Reader(myVcf):
				sType = siteType(record)
				if sType == 0:
					continue #skip over intronic sites
				pi = record.nucl_diversity
				piDict[sType].append(pi)
			outValues = [pac,sum(piDict["syn"]), sum(piDict["non"]), len(piDict["syn"]), len(piDict["non"])]
			out.write("	".join([str(x) for x in outValues])+"\n")

		elif myTest == "sfs":
			sfsDict = {"syn":[], "non":[]}
			for i in range(0,150): #make empty lists for the sfs
				sfsDict['syn'].append(0)
				sfsDict['non'].append(0)
			for record in vcf.Reader(myVcf):
				sType = siteType(record)
				if sType == 0:
					continue #skip over intronic sites
				myMaf = maf_downsample(record)	
				if myMaf == "NA":
					print(record)
					continue
				else:
					print(myMaf)
					sfsDict[sType][myMaf] += 1
			print(sfsDict)	



def maf(record): # calculates the maf
	af = record.aaf[0]
	if af > 0.5:
		maf = 1-af
	else:
		maf = af
	return(maf)

def maf_downsample(record): #calculates maf of a sample where n=160
	#remove any missing data
	noNas = []
	for call in record.samples:
		#print(call)
		#print(call.gt_type)
		if call.gt_type in [0,1,2]:
			noNas.append(call)
	if len(noNas) < 150:
		return("NA")  #not enough data
	mySample = random.sample(noNas,150) #randomly sample 150
	# how many are the rare homozygote?
	altCount = 0
	for call in mySample:
		altCount += call.gt_type

	if altCount > 149:
		maf = 300-altCount
	else:
		maf = altCount
	return(maf)

def siteType(record):#is it syn or nonsyn?
	siteType = 0
	for thing in record.INFO['EFF']:
		if thing.split('(')[0] == "NON_SYNONYMOUS_CODING":
			siteType = "syn"
		elif thing.split('(')[0] == "SYNONYMOUS_CODING":
			siteType = "non"
	return(siteType)
		

if __name__ == "__main__":
        __main__()





