gff = open('/cap1/emily.josephs/eQTL/data/gff/CR.gff','r')

out = open('../data/gff_gene_flatfile','w')

for line in gff:
	ent = line.split()
	scaffold, type, start, end, info = ent[0], ent[2], ent[3], ent[4], ent[8].split(';')
	if type != "mRNA":
		continue
	pac = info[0].split(':')[1]

	out.write('	'.join([scaffold, pac, start, end])+"\n")

	


