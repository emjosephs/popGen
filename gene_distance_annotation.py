# adds distance from nearest gene to Robert's annotation

import sys
if len(sys.argv) < 3:
	print('python gene_distance_annotation.py [all/any] [gff file] [summary input]')
	sys.exit()

def __main__():
	#make dict to hold distance values
	upstr = {}
	dwstr = {}
	for i in range(1,9):
		if sys.argv[1] == "any":
			upstr["scaffold_"+str(i)]=[5000]*19634518
			dwstr["scaffold_"+str(i)]=[5000]*19634518
		elif sys.argv[1] == "all":
			upstr["scaffold_"+str(i)]=[]
			dwstr["scaffold_"+str(i)]=[]
			for j in range(0,19634518):	
				upstr["scaffold_"+str(i)].append([])
				dwstr["scaffold_"+str(i)].append([])

	#read through gff
	gff = open(sys.argv[2],'r')
	for line in gff:
		if line[0] == "#":
			continue #header line
		lineDict = gff_parse(line)
		if lineDict['type'] != 'mRNA':
			continue
		if int(lineDict['scaf'].split('_')[1]) > 8:
			continue
		for i in range(0,5000):
			if sys.argv[1] == "any": #prints out distance from any tss or tes
				if int(lineDict['start']) > i: # check that start > i so start-i > 0, or else skip
					if int(upstr[lineDict['scaf']][int(lineDict['start'])-i]) > i: 
						upstr[lineDict['scaf']][int(lineDict['start'])-i] = i  # output list of distances from start
				if int(dwstr[lineDict['scaf']][int(lineDict['end'])+i]) > i: 
					dwstr[lineDict['scaf']][int(lineDict['end'])+i] = i  # output list of distances from iend
			

			elif sys.argv[1] == "all": #prints out distance from all tss/tes w/in 5kb	
				if int(lineDict['start']) > i: # this is for the gene clsoer to start of scaffold, check that start > i so start-i > 0, or else skip
					if lineDict['dir'] == '+':
						upstr[lineDict['scaf']][int(lineDict['start'])-i].append((lineDict['attrs']['pacid']+":"+str(i)))  #add gene name and distance upstream from tss
					else: #gene is reversed, site is downstream
						dwstr[lineDict['scaf']][int(lineDict['start'])-i].append((lineDict['attrs']['pacid']+":"+str(i)))  #add gene name and distance upstream from tss
				if lineDict['dir'] == '+':
					dwstr[lineDict['scaf']][int(lineDict['end'])+i].append((lineDict['attrs']['pacid']+":"+str(i)))  #add gene name and distance downstream from tss
				else:
					upstr[lineDict['scaf']][int(lineDict['end'])+i].append((lineDict['attrs']['pacid']+":"+str(i)))  #add gene name and distance downstream from tss
	
	gff.close()

	#read through annotatins and write new one
	#annot = open('/cap1/emily.josephs/genomes/scaf1_8.masked.new.downsampled.NCNC.summary','r')
	annot = open(sys.argv[3],'r')
	out = open(sys.argv[3]+"."+sys.argv[1]+'Distance','w')

	for line in annot:
		if line[0] == "#":
			out.write(line)
			continue
			
		annotNames = "CHROM  POS     REF     ALT     REF_NUMBER      ALT_NUMBER      TOTAL   SITE_TYPE       DIVERGENCE".split()
		lineDict = dict( zip( annotNames, line.split()[0:11] ) )

		out.write("	".join(line.split()[0:11]))
		if sys.argv[1] == "any":
			out.write( "	"+str( upstr[lineDict["CHROM"]][int(lineDict["POS"])] ) )
			out.write( "	"+str( dwstr[lineDict["CHROM"]][int(lineDict["POS"])] ) )
		elif sys.argv[1] == "all":
			for i in [upstr[lineDict["CHROM"]][int(lineDict["POS"])],dwstr[lineDict["CHROM"]][int(lineDict["POS"])]]:
				if len(i) > 0:
					out.write("	"+",".join(i))
				else:
					out.write("	NA")
		out.write("\n")

	annot.close()





def gff_parse(gffLine):
	#scaffold_1      JGI_gene        mRNA    750912  765975  .       -       .       ID=PAC:20892461;Name=Carubv10008059m;pacid=20892461;Parent=Carubv10008059m.g;
	names = "scaf	JGI_gene type start end . dir . attr".split()
	gffdict = dict( zip( names,gffLine.split()) )
	try:
		attrNames = [x.split("=")[0] for x in gffdict['attr'].split(';')[:-1]]
	except:
		print(gffLine)
		sys.exit()
	attrVals = [x.split("=")[1] for x in gffdict['attr'].split(';')[:-1]]
	gffdict['attrs'] = dict( zip( attrNames, attrVals ) )
	return(gffdict)

if __name__ == "__main__":
        
	#gff = open('/cap1/emily.josephs/eQTL/data/gff/CR.gff','r')
	#for line in gff:
	#	print(gff_parse(line))
	__main__()





