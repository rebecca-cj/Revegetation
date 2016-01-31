## Calculate Expected Heterozygosity by Chromosome for 11 chromsomes (based on Eucalypt genome) ##
# Input: individual population .frq files (vcftools output) ** N.B. assumes chromsome name is "scaffold_1" etc
# Output: 'HetByChr.csv' = one population per line with each line being Population, Chr Number, Mean Expected Heterozygosity, Standard Deviation, Number of loci
# Usage: python HeByChr.py

import glob
import numpy

Het={}
for fileX in glob.glob('*.frq'):
	tmpfile=open(fileX).readlines()
	sitename=fileX[:-4]
	Het[sitename]={}
	for chr in range(1,12):
		scaf="scaffold_"+str(chr)
   		tmphet=[]
    		for lineX in tmpfile:
       			raw=lineX.strip().split()
        		if raw[0]==scaf:
            			ref=raw[4].split(':')
            			alt=raw[5].split(':')
            			het=2*float(ref[1])*float(alt[1])
            			tmphet.append(het)
    		Het[sitename][scaf]=str((numpy.asarray(tmphet)).mean()) + "," + str((numpy.asarray(tmphet)).std(ddof=1)) + "," + str(len(tmphet))
    
outfile=open('HetByChr.csv','a')    
for sites in Het.keys():
	for chr in Het[sites].keys():
		outfile.write(sites + "," + chr[9:] + "," + Het[sites][chr] + "\n")	
outfile.close()
