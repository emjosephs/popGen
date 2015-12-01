#takes a list of arabidopsis genes, returns a list of capsella orthologs
import sys

def main():

	# read through ortho file and make a dict with at as keys and cg as values
	orthoDict = {}
	orthoFile = open(_args.ortho_file,'r')
	orthoFile.readline()
	for line in orthoFile:
		atName, cgName = line.split()[_args.ortho_at-1], line.split()[_args.ortho_cg - 1]
		orthoDict[atName] = cgName
	orthoFile.close()

	#read through list of AT genes and write out the CG names (with the AT names too?)
	atFile = open(_args.at_genes,'r')
	outFile = open(_args.output,'w')
	for line 




import argparse
_args = None
def parseArgs():
	parser = argparse.ArgumentParser(description = "AT to CG orthologs")
	parser.add_argument("-a", "--at_genes", type=str, help = "list of arabidopsis genes")
	parser.add_argument("-o", "--output", type=str, help = "output file")
	parser.add_argument("-c", "--column", type=int, default=1, help = "column with gene names, default is 1")
	parser.add_argument("-i", "--ortho_file", help = "list of orthologs, default is ~/genomes/CG_AT_PAC ", default = "/cap1/emily.josephs/genomes/CG_AT_PAC", type=str)
	parser.add_argument("-oa", "--ortho_at", help = "column in ortholog file with arabidopsis genes, default is 3", type=int, default=3)
	parser.add_argument("-oc", "--ortho_cg", help = "column in ortholog file with capsella genes, default is 4", type=int, default=4)
	

        global _args
        _args = parser.parse_args()
        sys.stderr.write(str(_args)+"\n")


if __name__ == "__main__":
        parseArgs()
        main()












