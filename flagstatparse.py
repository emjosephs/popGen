import os
danList = open('/cap1/emily.josephs/trans/bin/lists/dan_wgcna_list','r')
#for file in os.listdir('/cap2/GWAS_raw_data/bams/RNA_thirdflowcell/flagstat_out/'):
for line in danList:
       # indName = file.split('.')[0]
        myFile = open('/cap2/GWAS_raw_data/bams/RNA_thirdflowcell/flagstat_out/'+line.rstrip()+".flagstat",'r')
        allLines = myFile.readlines()
        total = allLines[0].split()[0]
        mapped = allLines[2].split()[0]

	print("    ".join([indName, total, mapped]))
