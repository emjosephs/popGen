import sys
import vcf

vcfList = open(sys.argv[1],'r')

for line in vcfList:
	myVcf = open(line.rstrip(),'r')
	
	for record in vcf.Reader(myVcf):





