import sys
import vcf

if len(sys.argv) < 4:
	print('python popgen_master.py [list of vcfs] [output file] [pi/sfs]')
	sys.exit()

def __main__():

	vcfList = open(sys.argv[1],'r')
	out = open(sys.argv[2],'w')
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
		out.write("pac	piS	piN	Ssite#	Nsite#\n")
		out.write("	".join([str(x) for x in outValues])+"\n")

						
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





