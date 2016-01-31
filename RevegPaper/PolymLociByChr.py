## Calculate percentage polymorphic loci by chromosome, for 11 chromosomes (based on Eucalypt genome) from population allele frequencies ##
# Input: individual population .frq files (vcftools output) ** N.B. assumes chromsome names in .frq are "scaffold_1" etc
# Output: PolymLociByChr.csv file = one population per line with each line being Site, Chr Number, Total number of loci, Number monomorphic loci, Number polymorphic loci, Percentage Polymorphic Loci
# Usage:  python PolymLociByChr.py

import glob

PLoci={}
for fileX in glob.glob('*.frq'):
	tmpfile=open(fileX).readlines()
	sitename=fileX[:-4]
	PLoci[sitename]={}
	for chr in range(1,12):
		scaf="scaffold_"+str(chr)
		total=0
		poly=0
		mono=0
		for lineX in tmpfile:
			raw=lineX.strip().split()
			if raw[0]==scaf:
				total=total+1
				ref=(raw[4].split(':'))[1]
				alt=(raw[5].split(':'))[1]
				if float(ref)==1 or float(alt)==1:
					mono=mono+1
				else:
					poly=poly+1
		PLoci[sitename][scaf]=str(total) + "," + str(mono) + "," + str(poly) + "," + str(float(poly)/float(total)*100)
		
outfile=open('PolymLociByChr.csv','a')
for sites in PLoci.keys():
	for chr in PLoci[sites].keys():
		outfile.write(sites + "," + chr[9:] + "," + PLoci[sites][chr] + "\n")
outfile.close()

