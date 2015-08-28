import os

for file in os.listdir('/cap2/GWAS_raw_data/bams/RNA_thirdflowcell/flagstat_out/'):
        indName = file.split('.')[0]
        myFile = open('/cap2/GWAS_raw_data/bams/RNA_thirdflowcell/flagstat_out/'+file,'r')
        allLines = myFile.readlines()
        total = allLines[0].split()[0]
        mapped = allLines[2].split()[0]

	print("    ".join([indName, total, mapped]))
