## Calculate Expected Heterozygosity by Chromosome for 11 chromsomes (based on Eucalypt genome) ##
#
#  Copyright 2016 Rebecca Jordan
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#####
# Input: individual population .frq files (vcftools output) ** N.B. assumes chromsome name is "scaffold_1" etc
# Output: 'HetByChr.csv' = one population per line with each line being Population, Chr Number, Mean Expected Heterozygosity, Standard Deviation, Number of loci
# Usage: python HeByChr.py
#####

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
